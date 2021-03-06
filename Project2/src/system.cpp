#include "system.h"
#include <cassert>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <algorithm>    // std::max
#include "sampler.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "Optimizers/optimizer.h"
#include "Math/random.h"

using namespace std;
using namespace arma;

/* Algorithm for standard Metropolis (brute force): change one coordinate at a time per MC cycle */
bool System::bruteForce() {
    int coor;
    double randu, s;
    arma::vec Xtrial = m_waveFunction->get_X();

    if(m_MCstep == 0){
        m_wfValue = m_waveFunction->evaluate(Xtrial);
        //        cout << "m_wfValue0 = " << m_wfValue << endl;
    }
    double currentWaveFunction = m_wfValue;
    // draw a random visible node
    coor = Random::nextInt(m_numberVisibleNodes);
    // compute weight
    s = Random::nextDouble();
    // compute change in coordinates from uniform distribution
    randu = getUniform(-1.0, 1.0);
    Xtrial[coor] += randu*m_stepLength;
    // evaluate wavefunction after change of coordinates
    double trialWaveFunction = m_waveFunction->evaluate(Xtrial);
    double ratio = (trialWaveFunction*trialWaveFunction)/(currentWaveFunction*currentWaveFunction); //Fix: can simplify expression to save cpu cycles

    // if true, allow the new state with adjusted positions
    if(s<ratio){
        m_waveFunction->set_X(Xtrial);
        m_wfValue = trialWaveFunction;

        return true;
    }

    // if false, reject the new state and change of coordinates
    else{
        m_wfValue = currentWaveFunction;

        return false;
    }
}


/* Algorithm for Metropolis-Hastings (importance sampling): change one coordinate at a time per MC cycle */
bool System::importanceSampling(){

    int coor;
    double randg, s, deriv;
    double QForceCurrent, QForceTrial;
    double D = 0.5;
    double GreensFunction = 0.0;
    double timestep = m_timesteps[0];
    double currentWaveFunction, trialWaveFunction;

    arma::vec Xcurrent = m_waveFunction->get_X();
    arma::vec Xtrial = m_waveFunction->get_X();
    if(m_MCstep == 0){
        m_wfValue = m_waveFunction->evaluate(Xcurrent);
    }
    currentWaveFunction = m_wfValue;

    // draw a random visible node
    coor = Random::nextInt(m_numberVisibleNodes);
    // compute weight
    s = Random::nextDouble();

    // compute current quantum force and wavefunction
    deriv = m_waveFunction->computeDerivative_analytic(Xcurrent, coor);
    QForceCurrent = 2*deriv;

    // compute change in coordinates from uniform distribution
    randg = getGaussian(0, 1);
    Xtrial[coor] += D*QForceCurrent*timestep + randg*sqrt(timestep);

    // compute trial quantum force and wavefunction after change of coordinates
    trialWaveFunction = m_waveFunction->evaluate(Xtrial);
    deriv = m_waveFunction->computeDerivative_analytic(Xtrial, coor);
    QForceTrial = 2*deriv;

    //Greens ratio
    double part1 = Xcurrent[coor] - Xtrial[coor] - timestep*QForceTrial;
    double part2 = Xtrial[coor] - Xcurrent[coor] - timestep*QForceCurrent;
    GreensFunction = exp(-(part1*part1 - part2*part2)/(4*D*timestep));

    double probRatio = GreensFunction*trialWaveFunction*trialWaveFunction/(currentWaveFunction*currentWaveFunction);

    // if true, allow the new state with adjusted positions
    if(s<probRatio){
        m_waveFunction->set_X(Xtrial);
        m_wfValue = trialWaveFunction;

        return true;
    }

    // if false, reject the new state and change of coordinates
    else{
        m_wfValue = currentWaveFunction;

        return false;
    }
}

/* Algorithm for Gibbs sampling */
bool System::Gibbs(){
    double randu, P, sigma;
    arma::vec O;

    sigma = m_waveFunction->getSigma();
    arma::vec X = m_waveFunction->get_X();
    arma::vec m_h = m_waveFunction->get_h();
    arma::vec m_a = m_waveFunction->get_a();
    arma::vec m_b = m_waveFunction->get_b();
    arma::mat m_w = m_waveFunction->get_w();

    // Get probability, P(h=1|x), of hidden values to equal 1, given the logistic sigmoid function
    // Set hidden values equal to 0 if probability less than a randuom uniform variable
    O = m_b + ((X.t()*m_w).t()) / (sigma*sigma);
    for (int j=0; j < m_numberHiddenNodes; j++){
        randu = getUniform(0, 1);
        P = 1.0/(1+exp(-O[j]));     // probability from logistic sigmoid
        if (P < randu){ m_h[j] = 0; }
        else{ m_h[j] = 1; }
    }

    // Set the new positions according to the hidden nodes
    double randn, x_mean;
    arma::vec new_pos, wh;
    new_pos.zeros(m_numberVisibleNodes);

    wh = m_w*m_h;
    for (int i=0; i<m_numberVisibleNodes; i++){
        x_mean = m_a[i] + wh[i];
        new_pos[i] = getGaussian(x_mean, sigma);
    }

//    cout << "X before: " << endl; m_waveFunction->get_X().print();
    m_waveFunction->set_X(new_pos);
//    cout << "X after: " << endl; m_waveFunction->get_X().print();

    return true;
}


void System::runMetropolisSteps(int RBM_cycles, std::vector<int> numberOfMetropolisSteps) {

    m_timestep = 0;
    m_numberOfMetropolisSteps = numberOfMetropolisSteps[0];
    m_RBMcycles = RBM_cycles;
    m_MCstep = 0;
    m_acceptedSteps = 0.0;
    m_sampler = new Sampler(this);
    m_sampler->setNumberOfMetropolisSteps(m_numberOfMetropolisSteps);
    double m_initialization = m_waveFunction->getInitializationInterval();

    std::random_device seed;
    mt19937_64 gen(seed());
    m_randomengine = gen;

    // sampling cycles
    for (int cycle = 0; cycle < m_RBMcycles; cycle++){
        cout << "---------------------------------------" << endl;
        cout << "RBM cycle = " << cycle << endl;
        cout << "---------------------------------------" << endl;
        m_MCstep = 0;
        m_sampleStep = 0;
        // run Metropolis algorithm
        m_acceptedSteps = 0;
        for (int i=0; i < m_numberOfMetropolisSteps; i++) {

            bool acceptedStep;

            // set the solver
            if (m_solver=="bruteForce"){ acceptedStep = bruteForce(); }
            else if (m_solver=="importance"){ acceptedStep = importanceSampling(); }
            else if (m_solver=="gibbs"){ acceptedStep = Gibbs(); }


            // sample steps
            if (acceptedStep == true){
                m_acceptedSteps++;
            }

            // allow for equilibration of energy (~10% of Metropolis steps)
            if (i >= m_equilibrationFraction){        // m_equilibrationFraction - 100
                m_sampler->sample(acceptedStep);

                // Only interested in sampling the final optimisation cycle
                if (cycle == RBM_cycles - 1 ){
                    string filename_blocking = "../data/g/blocking/2g-Initialization-";
                    filename_blocking.append(to_string(m_initialization));
                    filename_blocking.append("_blockingSteps_");
                    m_sampler->writeStepToFile(m_sampleStep, m_sampleStep, filename_blocking);
                }
                m_sampleStep++;
            }
            m_MCstep++;
        }
        // write energies of RBM cycles to file
        string filename_RBM = "../data/g/RBM/2g-Initialization-";
        filename_RBM.append(to_string(m_initialization));
        filename_RBM.append("_RBMcycles_");
        m_sampler->writeToFile(cycle, cycle, filename_RBM);
        m_sampler->computeAverages();
        m_sampler->printOutputToTerminal();
        m_acceptedSteps_ratio = m_acceptedSteps/((double) m_numberOfMetropolisSteps);
        cout << "Acceptance rate: " << m_acceptedSteps_ratio << "\n" << endl;

        // optimize weights at end of cycle
        m_sampler->optimizeWeights();

        m_RBMstep++;
    }
}


void System::setNumberParticles(int nPart) {
    m_numberParticles = nPart;
}

void System::setNumberDimensions(int nDim) {
    m_numberDimensions = nDim;
}

void System::setStepLength(double stepLength) {
    assert(stepLength >= 0);
    m_stepLength = stepLength;
}

void System::setEquilibrationFraction(double equilibrationFraction) {
    assert(equilibrationFraction >= 0);
    m_equilibrationFraction = equilibrationFraction;
}

void System::setHamiltonian(Hamiltonian* hamiltonian) {
    m_hamiltonian = hamiltonian;
}

void System::setWaveFunction(WaveFunction* waveFunction) {
    m_waveFunction = waveFunction;
}

void System::setOptimizer(Optimizer* optimizer){
    m_optimizer = optimizer;
}
