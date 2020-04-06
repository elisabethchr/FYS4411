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

bool System::metropolisStep() {
    /* Draw a particle at random and change its position by a random amount.
     * Evaluate the wavefunctoin at this new position and check if step to the new state
     * is accepted.
     */

    if(m_MCstep == 0){
        m_wfValue = m_waveFunction->evaluate(m_particles);
        cout << "m_wfValue = " << m_wfValue << endl;
    }

    int random_i;
    double random_d, change, s;

    // draw random particle i
    random_i = Random::nextInt(m_numberOfParticles);
    // compute weight
    s = Random::nextDouble();

    arma::vec change_vec(m_numberOfDimensions);
    arma::vec rand(m_numberOfDimensions);

    change_vec.zeros();
    //    std::vector<double> alphas =  m_waveFunction->getParameters();

    // evaluate wavefunction at old position
    double oldWaveFunction = m_waveFunction->evaluate(m_particles);
    //double oldWaveFunction = m_wfValue;

    // change positions
    for(int dim=0; dim<m_numberOfDimensions; dim++){
        random_d = getUniform(-1.0, 1.0);
        rand[dim] = random_d;
        change = m_stepLength*random_d;       // between -0.5 and 0.5?

        change_vec[dim] = change;
        m_particles[random_i]->adjustPosition(change, dim);
    }

    // evaluate wavefunction at new position
    double newWaveFunction = m_waveFunction->evaluate(m_particles);

    double ratio = (newWaveFunction*newWaveFunction)/(oldWaveFunction*oldWaveFunction); //Fix: can simplify expression to save cpu cycles

    // if true, allow the new state with adjusted positions
    if(s<=ratio){
        m_wfValue = newWaveFunction;
        return true;
    }
    // if false, reject the new state and reset positions
    else{
        m_wfValue = oldWaveFunction;
        for(int dim=0; dim<m_numberOfDimensions; dim++){
            m_particles[random_i]->adjustPosition(-change_vec[dim], dim);
        }
        return false;
    }
}


bool System::importanceSampling(){
    /* Draw a particle at random, change its position according to a normal distribution and calculate the new wavefunction.
     * Calculate the corresponding quantum forces before and after the move, and compute the Green's function G as a function of the positions and
     * the quantum forces. The ratio is then r = G*wfNew^2/wfOld^2.
    */

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
    QForceOld = m_hamiltonian->computeQuantumForce(m_particles, random_i); // /oldWaveFunction; // (dividing by m_wfValue in computeQuantumForce for harmonic oscillator)

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
    QForceNew = m_hamiltonian->computeQuantumForce(m_particles, random_i); // /newWaveFunction; // (dividing by m_wfValue in computeQuantumForce for harmonic oscillator)

    // Compute Green's function by looping over all dimensions, where m_stepLength ~= timestep
    for (int j=0; j<m_numberOfDimensions; j++){
        GreensFunction += 0.5 * (QForceOld[j] + QForceNew[j]) *(D*timestep*0.5 * (QForceOld[j] - QForceNew[j]) - posNew[random_i]->getPosition()[j] + posOld[random_i]->getPosition()[j]);
    }

    GreensFunction = exp(GreensFunction);

    double ratio = GreensFunction*newWaveFunction*newWaveFunction/(oldWaveFunction*oldWaveFunction); //Fix: can simplify expression to save cpu cycles

    if(s<=ratio){
        //m_wfValue = newWaveFunction;
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
}


void System::runMetropolisSteps(std::vector<int> numberOfMetropolisSteps) {
    m_alpha = 0;
    for (int alpha=0; alpha<m_waveFunction->getParameters().size(); alpha++){
        cout << "\n m_alpha = " << m_alpha << ", " << "alpha = " << m_waveFunction->getParameters()[m_alpha] << endl;
//        cout << "m_beta = " << m_waveFunction->getParametersBeta()[0] << endl;
        cout << "m_gamma = " << m_waveFunction->getGamma() << endl;
        cout << "m_omega = " << m_hamiltonian->getOmega() << endl;

        //        for(int MC=0; MC<numberOfMetropolisSteps.size(); MC++){
        m_timestep = 0;
        m_numberOfMetropolisSteps   = numberOfMetropolisSteps[0];
        cout << "\n m_numberOfMetropolisStep: " << "2^" << log2(m_numberOfMetropolisSteps) << " = " << m_numberOfMetropolisSteps << endl;
        for (int dt=0; dt<m_timesteps.size(); dt++){
            cout << "\n m_timestep = " << m_timestep << ", " << "timestep = " << m_timesteps[0] << endl;
            m_stepMetropolis = 0.0;
            m_stepImportance = 0.0;
            m_acceptedSteps = 0.0;
            m_MCstep = 0;
            m_particles                 = m_initialState->getParticles();
            m_sampler                   = new Sampler(this); //Remove later: (this) points to the system object from which  "runMetropolisSteps" is called.
            //m_numberOfMetropolisSteps   = numberOfMetropolisSteps[MC];
            m_sampler->setNumberOfMetropolisSteps(m_numberOfMetropolisSteps);
            int steps = 0;

            std::random_device i;
            mt19937_64 gen(i());
            m_seed = gen;

            for (int i=0; i < m_numberOfMetropolisSteps; i++) {
                bool acceptedStep;

                // set the solver
                if(m_solver==true){ acceptedStep = metropolisStep(); }
                else if (m_solver==false){ acceptedStep = importanceSampling(); }

                if(acceptedStep == true){
                    m_acceptedSteps++;
                    if(i >= m_equilibrationFraction - 100){
                        m_sampler->sample(acceptedStep);

                        steps++;
                    }
                }
//                m_sampler->computeAverages();

                //if (i>= m_equilibrationFraction){ m_sampler->writeStepToFile(i, i); }

                m_stepMetropolis++;
                m_MCstep++;
            }
            m_sampler->computeAverages();
            m_sampler->printOutputToTerminal();
           // m_sampler->writeAlphaToFile();
            m_acceptedSteps_ratio = m_acceptedSteps/((double) m_numberOfMetropolisSteps);
            cout << "Acceptance rate: " << m_acceptedSteps_ratio << endl;

            m_energy = m_sampler->getEnergy();
            m_derivativeE = m_sampler->getEnergyDerivative();

            m_timestep++;
        }
        //}
        m_alpha++;
    }
}



void System::setNumberOfParticles(int numberOfParticles) {
    m_numberOfParticles = numberOfParticles;
}

void System::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
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
