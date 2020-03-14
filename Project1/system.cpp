#include "system.h"
#include <cassert>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"

using namespace std;
using namespace arma;

bool System::metropolisStep() {
    /* Perform the actual Metropolis step: Choose a particle at random and
     * change its position by a random amount, and check if the step is
     * accepted by the Metropolis test (compare the wave function evaluated
     * at this new position with the one at the old position).
     */
    /*
     * Look at/draw one particle at a time. Allocate matrices which contain the position of the particle.
    */

    double random_d; //Random dimension
    int random_i;
    double s;
    random_i = Random::nextInt(m_numberOfParticles);
    s = Random::nextDouble();
    double dx;

    arma::vec dx_vec(m_numberOfDimensions); dx_vec.zeros();
    std::vector<double> alphas =  m_waveFunction->getParameters();

//    double dx = m_stepLength*(random_d - 0.5);

    double oldWaveFunction = m_waveFunction->evaluate(m_particles);

    for(int dim=0; dim<m_numberOfDimensions; dim++){
        random_d = Random::nextDouble();
        dx = m_stepLength*(random_d - 0.5);
        dx_vec[dim] = dx;
        m_particles[random_i]->adjustPosition(dx, dim);
    }

    double newWaveFunction = m_waveFunction->evaluate(m_particles);

    double ratio = newWaveFunction*newWaveFunction/(oldWaveFunction*oldWaveFunction); //Fix: can simplify expression to save cpu cycles
    if(s<=ratio){
        return true;
    }else{    for(int dim=0; dim<m_numberOfDimensions; dim++){
            m_particles[random_i]->adjustPosition(-dx_vec[dim], dim);
        }
        return false;
    }
}

bool System::importanceSampling(){
    double random_d; //Random dimension
    int random_i;
    double s;
    random_i = Random::nextInt(m_numberOfParticles);
    s = Random::nextDouble();
    double dx;

    arma::vec dx_vec(m_numberOfDimensions); dx_vec.zeros();
    std::vector<double> alphas =  m_waveFunction->getParameters();

//    double dx = m_stepLength*(random_d - 0.5);

    double oldWaveFunction = m_waveFunction->evaluate(m_particles);

    for(int dim=0; dim<m_numberOfDimensions; dim++){
        random_d = Random::nextDouble();
        dx = m_stepLength*(random_d - 0.5);
        dx_vec[dim] = dx;
        m_particles[random_i]->adjustPosition(dx, dim);
    }

    double newWaveFunction = m_waveFunction->evaluate(m_particles);

    double ratio = newWaveFunction*newWaveFunction/(oldWaveFunction*oldWaveFunction); //Fix: can simplify expression to save cpu cycles
    if(s<=ratio){
        return true;
    }else{    for(int dim=0; dim<m_numberOfDimensions; dim++){
            m_particles[random_i]->adjustPosition(-dx_vec[dim], dim);
        }
        return false;
    }
}

void System::runMetropolisSteps(int numberOfMetropolisSteps) {
    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this); //Remove later: (this) points to the system object from which  "runMetropolisSteps" is called.
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);
    double steps = 0;

    for (int i=0; i < numberOfMetropolisSteps; i++) {
        for (int j=0; j < getNumberOfParticles(); j++){

            /*Update positions of particles*/
            bool acceptedStep = metropolisStep();
//            bool acceptedStep = importanceSampling();
            if(acceptedStep == true){

                m_sampler->sample(acceptedStep);
//                if ((i%100==0) && (i != 0)){steps ++;}
                /* Here you should sample the energy (and maybe other things using
         * the m_sampler instance of the Sampler class. Make sure, though,
         * to only begin sampling after you have let the system equilibrate
         * for a while. You may handle this using the fraction of steps which
         * are equilibration steps; m_equilibrationFraction.
         */

            }
        }
        // write only every 10 value
        if ((m_numberOfMetropolisSteps >= 1e4) && (i%10==0) && (i != 0)){
            m_sampler->writeToFile(i, steps);
            steps++;
        }
        if ((m_numberOfMetropolisSteps < 1e4) && (i%1==0) && (i != 0)){
            m_sampler->writeToFile(i, steps);
            steps++;
        }
        //m_sampler->writeToFile(i);
    }
    m_sampler->computeAverages();
    m_sampler->printOutputToTerminal();
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


