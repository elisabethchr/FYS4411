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
    /* Perform the actual Metropolis step: Choose a particle at random and
     * change its position by a random amount, and check if the step is
     * accepted by the Metropolis test (compare the wave function evaluated
     * at this new position with the one at the old position).
     */
    /*
     * Look at/draw one particle at a time. Allocate matrices which contain the position of the particle.
    */

    if(m_MCstep == 0){
    cout << "evaluating wfValue first MC_step" << endl;
    m_wfValue = m_waveFunction->evaluate(m_particles);
    cout << "First step: " << m_wfValue << endl;
    }

    double random_d; //Random dimension
    int random_i;
    double s;
    random_i = Random::nextInt(m_numberOfParticles);
    s = Random::nextDouble();
    double dx;

    arma::vec dx_vec(m_numberOfDimensions);
    arma::vec rand(m_numberOfDimensions);

    dx_vec.zeros();
    std::vector<double> alphas =  m_waveFunction->getParameters();

    // evaluate wavefunction before adjusting positions
    double oldWaveFunction = m_waveFunction->evaluate(m_particles);
//    double oldWaveFunction = m_wfValue;

//    if(oldWaveFunction != oldWaveFunction_test){
//    cout << "oldWavefunction = " << setprecision(8) << oldWaveFunction << ", m_wfValue = " << setprecision(8) << oldWaveFunction_test << endl;
//    }

    // testing
    for(int dim=0; dim<m_numberOfDimensions; dim++){
        random_d = getUniform(-1.0, 1.0);
        rand[dim] = random_d;
        dx = m_stepLength*random_d;       // between -0.5 and 0.5?
//        cout << "dx = " << dx << endl;
        dx_vec[dim] = dx;
        m_particles[random_i]->adjustPosition(dx, dim);
    }
/*
    // original
    for(int dim=0; dim<m_numberOfDimensions; dim++){
        random_d = Random::nextDouble();
//        cout << "random::nextDouble = " << random_d << endl;
        dx = m_stepLength*(random_d-0.5);       // between -0.5 and 0.5
        dx_vec[dim] = dx;
        m_particles[random_i]->adjustPosition(dx, dim);
    }
*/
    // evaluate wavefunction after adjusting positions
    double newWaveFunction = m_waveFunction->evaluate(m_particles);

    double ratio = (newWaveFunction*newWaveFunction)/(oldWaveFunction*oldWaveFunction); //Fix: can simplify expression to save cpu cycles
    // if true, allow the new state with adjusted positions
    if(s<=ratio){
        //        cout << "True: " << newWaveFunction << endl;
        m_wfValue = newWaveFunction;
        return true;
    }
    // if false, reject the new state and reset positions
    else{
//        cout << "False: " << oldWaveFunction << endl;
        m_wfValue = oldWaveFunction;
        for(int dim=0; dim<m_numberOfDimensions; dim++){
            m_particles[random_i]->adjustPosition(-dx_vec[dim], dim);
        }
        return false;
    }
}


bool System::importanceSampling(){
    int random_i;
    double s, random_d, dx;
    double D = 0.5;
    double GreensFunction = 0.0;
    double h = 1e-4;
    double timestep = 0.001;
    double oldWaveFunction, newWaveFunction, quantumForceOld, quantumForceNew;
    random_i = Random::nextInt(m_numberOfParticles);
    s = Random::nextDouble();

    std::vector<double> alphas =  m_waveFunction->getParameters();

    arma::mat dx_mat(m_numberOfParticles, m_numberOfDimensions); dx_mat.zeros();
    arma::mat QForceOld(m_numberOfParticles, m_numberOfDimensions); QForceOld.zeros();
    arma::mat QForceNew(m_numberOfParticles, m_numberOfDimensions); QForceNew.zeros();

    // Old position
    std::vector<Particle *> posOld = m_particles;
    oldWaveFunction = m_waveFunction->evaluate(m_particles);
    QForceOld = m_hamiltonian->computeQuantumForce(m_particles)/oldWaveFunction;

    // Move a random distance dx in all dimensions
    for (int i =0; i<m_numberOfParticles; i++){
        for (int dim=0; dim<m_numberOfDimensions; dim++){
            //            random_d = Random::nextDouble();
            //            dx = m_stepLength*(random_d - 0.5);
            //            dx_mat(i, dim) = dx;
            m_particles[i]->adjustPosition(dx_mat(i, dim), dim);
        }
    }

    // New position
    std::vector<Particle *> posNew = m_particles;
    newWaveFunction = m_waveFunction->evaluate(m_particles);
    QForceNew = m_hamiltonian->computeQuantumForce(m_particles)/newWaveFunction;


    // Compute Green's function by looping over all particles and dimensions, where m_stepLength ~= timestep
    for (int i=0; i<m_numberOfParticles; i++){
        for (int j=0; j<m_numberOfDimensions; j++){
            GreensFunction += 0.5*(QForceOld(i, j) + QForceNew(i, j))*(D*timestep*0.5*(QForceOld(i, j) - QForceNew(i, j)) - posNew[i]->getPosition()[j] + posOld[i]->getPosition()[j]);
        }
    }

    GreensFunction = exp(GreensFunction);

    double ratio = GreensFunction*newWaveFunction*newWaveFunction/(oldWaveFunction*oldWaveFunction); //Fix: can simplify expression to save cpu cycles

    if(s<=ratio){
        m_stepImportance++;
        return true;
    }
    else{
        for (int i=0; i<m_numberOfParticles; i++){
            for(int dim=0; dim<m_numberOfDimensions; dim++){
                m_particles[i]->adjustPosition(-dx_mat(i, dim), dim);
            }
        }
        return false;
    }
}


void System::runMetropolisSteps(int numberOfMetropolisSteps) {

    for (int alpha=0; alpha<m_waveFunction->getParameters().size(); alpha++){
        m_stepMetropolis = 0.0;
        m_stepImportance = 0.0;
        m_MCstep = 0;
        cout << "\n m_alpha = " << m_alpha << ", " << "alpha = " << m_waveFunction->getParameters()[m_alpha] << endl;
        m_particles                 = m_initialState->getParticles();
        m_sampler                   = new Sampler(this); //Remove later: (this) points to the system object from which  "runMetropolisSteps" is called.
        m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
        m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);
        bool type                   = getCalculation();

        unsigned __int64 i;
        i = __rdtsc();
        mt19937_64 gen(i);
        m_seed = gen;

        for (int i=0; i < numberOfMetropolisSteps; i++) {

            bool acceptedStep = metropolisStep();
            //            bool acceptedStep1 = importanceSampling();

            if(acceptedStep == true){
                m_stepMetropolis++;
                if(i > m_equilibrationFraction*numberOfMetropolisSteps){
                    m_sampler->sample(acceptedStep);
                }
                //                if ((i%100==0) && (i != 0)){steps ++;}
                //          Here you should sample the energy (and maybe other things using
                //          the m_sampler instance of the Sampler class. Make sure, though,
                //          to only begin sampling after you have let the system equilibrate
                //          for a while. You may handle this using the fraction of steps which
                //          are equilibration steps; m_equilibrationFraction.


            }
            //      }
            // write only every 10 value

            //        if ((m_numberOfMetropolisSteps >= 1e4) && (i%10==0) && (i != 0)){
            //            m_sampler->writeStepToFile(i, steps);
            //            steps++;
            //        }
            //        if ((m_numberOfMetropolisSteps < 1e4) && (i%1==0) && (i != 0)){
            //            m_sampler->writeStepToFile(i, steps);
            //            steps++;
            //        }

            //m_sampler->writeStepToFile(i);
        m_MCstep++;
        }
        m_sampler->computeAverages();
        m_sampler->printOutputToTerminal();

        m_sampler->writeAlphaToFile();

        cout << "Acceptance rate Metropolis: " << m_stepMetropolis/((double) m_numberOfMetropolisSteps) << endl;
        cout << "Acceptance rate importance sampling: " << m_stepImportance/((double) m_numberOfMetropolisSteps) << endl;

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


