#include "system.h"
#include <cassert>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <algorithm>    // std::max
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
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
//        cout << "m_wfValue = " << m_wfValue << endl;
    }

    // draw a random coordinate
    coor = Random::nextInt(m_numberVisibleNodes);        // number of visible nodes

    // compute weight
    s = Random::nextDouble();

    // compute change in coordinates from uniform distribution
    randu = getUniform(-1.0, 1.0);

//    cout << "Current Xtrial: " << endl;
//    for (int i=0; i<Xtrial.size(); i++){ cout << Xtrial[i] << endl; }

    // evaluate wavefunction before and after change of coordinates
    double oldWaveFunction = m_waveFunction->evaluate(Xtrial);
    Xtrial[coor] += randu*m_stepLength;
    double newWaveFunction = m_waveFunction->evaluate(Xtrial);

    double ratio = (newWaveFunction*newWaveFunction)/(oldWaveFunction*oldWaveFunction); //Fix: can simplify expression to save cpu cycles

//    cout << "New Xtrial: " << endl;
//    for (int i=0; i<Xtrial.size(); i++){ cout << Xtrial[i] << endl; }


    // if true, allow the new state with adjusted positions
    if(s<ratio){
        m_waveFunction->set_X(Xtrial);
        m_wfValue = newWaveFunction;

        return true;
    }

    // if false, reject the new state and reset positions
    else{
        Xtrial[coor] -= randu*m_stepLength;
        m_waveFunction->set_X(Xtrial);
        m_wfValue = oldWaveFunction;

        return false;
    }
}


/* Algorithm for Metropolis-Hastings (importance sampling): change one coordinate at a time per MC cycle */
bool System::importanceSampling(){
/*
    int random_i;
    double s, change;
    double D = 0.5;
    double GreensFunction = 0.0;
    double timestep = m_timesteps[0];
    double oldWaveFunction, newWaveFunction;
    arma::vec change_vec(m_numberOfDimensions); change_vec.zeros();
    arma::vec QForceOld(m_numberOfDimensions); QForceOld.zeros();
    arma::vec QForceNew(m_numberOfDimensions); QForceNew.zeros();

    // draw random particle i
    random_i = Random::nextInt(m_numberOfParticles);
    // compute weight
    s = Random::nextDouble();


    // Old position
    std::vector<Particle *> posOld = m_particles;
    oldWaveFunction = m_waveFunction->evaluate(m_particles);
    m_wfValue = oldWaveFunction;
    QForceOld = m_hamiltonian->computeQuantumForce(m_particles, random_i); // (dividing by m_wfValue in computeQuantumForce for harmonic oscillator)

    // Move a random distance in every dimension according to a Gaussian distribution
    for (int dim=0; dim<m_numberOfDimensions; dim++){
        change = getGaussian(0, 1)*pow(timestep, 0.5) + QForceOld[dim]*timestep*D;

        m_particles[random_i]->adjustPosition(change, dim);
        change_vec[dim] = change;
    }

    // New position
    std::vector<Particle *> posNew = m_particles;
    newWaveFunction = m_waveFunction->evaluate(m_particles);
    m_wfValue = newWaveFunction;
    QForceNew = m_hamiltonian->computeQuantumForce(m_particles, random_i); // (dividing by m_wfValue in computeQuantumForce for harmonic oscillator)

    // Compute Green's function by looping over all dimensions, where m_stepLength ~= timestep
    for (int j=0; j<m_numberOfDimensions; j++){
        GreensFunction += 0.5 * (QForceOld[j] + QForceNew[j]) *(D*timestep*0.5 * (QForceOld[j] - QForceNew[j]) - posNew[random_i]->getPosition()[j] + posOld[random_i]->getPosition()[j]);
    }

    GreensFunction = exp(GreensFunction);

    double ratio = GreensFunction*newWaveFunction*newWaveFunction/(oldWaveFunction*oldWaveFunction); //Fix: can simplify expression to save cpu cycles

    if(s<=ratio){
        m_stepImportance++;
        return true;
    }
    else{
        m_wfValue = oldWaveFunction;
        for(int dim=0; dim<m_numberOfDimensions; dim++){
            m_particles[random_i]->adjustPosition(-change_vec(dim), dim);
        }
        return false;
    }
*/
    return true;
}


void System::runMetropolisSteps(std::vector<int> numberOfMetropolisSteps) {

    m_timestep = 0;
    m_numberOfMetropolisSteps = numberOfMetropolisSteps[0];
    m_stepMetropolis = 0.0;
    m_stepImportance = 0.0;
    m_acceptedSteps = 0.0;
    m_MCstep = 0;
    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this);
    m_sampler->setNumberOfMetropolisSteps(m_numberOfMetropolisSteps);

    int steps = 0;

    std::random_device seed;
    mt19937_64 gen(seed());
    m_randomengine = gen;

    // run Metropolis algorithm
    for (int i=0; i < m_numberOfMetropolisSteps; i++) {

        bool acceptedStep;

        // set the solver
        if (m_solver==true){ acceptedStep = bruteForce(); }
        else if (m_solver==false){ acceptedStep = importanceSampling(); }

        // sample accepted steps
        if (acceptedStep == true){
            m_acceptedSteps++;
            // allow for equilibration of energy (~5% of Metropolis steps)
            if (i >= m_equilibrationFraction - 100){
                m_sampler->sample(acceptedStep);
                steps++;
                m_stepMetropolis++;
            }
        }

        if (i>= m_equilibrationFraction){ m_sampler->writeStepToFile(m_stepMetropolis, i); }

        // Only interested in sampling the final optimization cycle
        if (i == m_numberOfMetropolisSteps - 1){ m_sampler->writeToFile(); }

        m_MCstep++;
    }
    m_sampler->computeAverages();
    m_sampler->printOutputToTerminal();
//    m_sampler->writeAlphaToFile();
    m_acceptedSteps_ratio = m_acceptedSteps/((double) m_numberOfMetropolisSteps);
    cout << "Acceptance rate: " << m_acceptedSteps_ratio << endl;
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

void System::setInitialState(InitialState* initialState) {
    m_initialState = initialState;
}
