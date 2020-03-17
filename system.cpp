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
 //   cout << "evaluating wfValue first MC_step" << endl;
    m_wfValue = m_waveFunction->evaluate(m_particles);
 //   cout << "First step: " << m_wfValue << endl;
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
        m_wfValue = newWaveFunction;
        return true;
    }
    // if false, reject the new state and reset positions
    else{
        m_wfValue = oldWaveFunction;
        for(int dim=0; dim<m_numberOfDimensions; dim++){
            m_particles[random_i]->adjustPosition(-dx_vec[dim], dim);
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
    double s, random_d, change, gauss;
    double D = 0.5;
    double GreensFunction = 0.0;
    double h = 1e-4;
    double timestep = 0.01;
    double oldWaveFunction, newWaveFunction, quantumForceOld, quantumForceNew;
    random_i = Random::nextInt(m_numberOfParticles);
    s = Random::nextDouble();

    arma::vec change_vec(m_numberOfDimensions); change_vec.zeros();
    arma::vec QForceOld(m_numberOfDimensions); QForceOld.zeros();
    arma::vec QForceNew(m_numberOfDimensions); QForceNew.zeros();

    // Old position
    std::vector<Particle *> posOld = m_particles;
    oldWaveFunction = m_waveFunction->evaluate(m_particles);
//    cout << "oldWavefunction = " << oldWaveFunction << endl;
    QForceOld = m_hamiltonian->computeQuantumForce(m_particles, random_i)/oldWaveFunction;

//    cout << "QForceOld" << QForceOld << endl;

    // Move a random distance in every dimensions
//    for (int i =0; i<m_numberOfParticles; i++){
        for (int dim=0; dim<m_numberOfDimensions; dim++){
            gauss = getGaussian(0, 1);
            change = gauss*pow(timestep, 0.5)+QForceOld[dim]*timestep*D;
//            cout << "Gauss" << gauss << endl;
            m_particles[random_i]->adjustPosition(change, dim);
//            cout << "change = " << change << endl;
        }
//    }

    // New position
    std::vector<Particle *> posNew = m_particles;
    newWaveFunction = m_waveFunction->evaluate(m_particles);
    QForceNew = m_hamiltonian->computeQuantumForce(m_particles, random_i)/newWaveFunction;

//    cout << "newWavefunction = " << newWaveFunction << endl;

//    cout << "QForceNew" << QForceNew << endl;

    // Compute Green's function by looping over all dimensions, where m_stepLength ~= timestep
//    for (int i=0; i<m_numberOfParticles; i++){
        for (int j=0; j<m_numberOfDimensions; j++){
            GreensFunction += 0.5*(QForceOld[j] + QForceNew[j])*(D*timestep*0.5*(QForceOld[j] - QForceNew[j]) - posNew[random_i]->getPosition()[j] + posOld[random_i]->getPosition()[j]);
        }
//    }

    GreensFunction = exp(GreensFunction);
//    cout << GreensFunction << endl;

    double ratio = GreensFunction*newWaveFunction*newWaveFunction/(oldWaveFunction*oldWaveFunction); //Fix: can simplify expression to save cpu cycles
//cout << "ratio:     " <<ratio<<endl;
    if(s<=ratio){
        m_stepImportance++;
        m_wfValue = newWaveFunction;
        return true;
    }
    else{
//        for (int i=0; i<m_numberOfParticles; i++){
            for(int dim=0; dim<m_numberOfDimensions; dim++){
                m_particles[random_i]->adjustPosition(-change_vec(dim), dim);
            }
//        }
        m_wfValue = oldWaveFunction;
        return false;
    }
}


void System::runMetropolisSteps(std::vector<int> numberOfMetropolisSteps, bool bruteForce) {


    for (int alpha=0; alpha<m_waveFunction->getParameters().size(); alpha++){
           cout << "\n m_alpha = " << m_alpha << ", " << "alpha = " << m_waveFunction->getParameters()[m_alpha] << endl;
           for (int dt=0; dt<m_timesteps.size(); dt++){
               cout << "\n m_timestep = " << m_timestep << ", " << "timestep = " << m_timesteps[m_timestep] << endl;

               m_nStepIndex = 0;
               for (int n=0; n < numberOfMetropolisSteps.size(); n++){
                   cout << n << endl;
               m_stepMetropolis = 0.0;
               m_stepImportance = 0.0;
               m_acceptedSteps = 0.0;
               m_MCstep = 0;
               m_particles                 = m_initialState->getParticles();
               m_sampler                   = new Sampler(this); //Remove later: (this) points to the system object from which  "runMetropolisSteps" is called.
//               m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
               m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);
               int steps = 0;

               std::random_device i;
               mt19937_64 gen(i());
               m_seed = gen;
               cout << "numberOfSteps size "<<numberOfMetropolisSteps.size()<<endl;


               for (int i=0; i < numberOfMetropolisSteps[n]; i++) {
                   bool acceptedStep;

                   if(bruteForce){
                      acceptedStep = metropolisStep();
                   }else{
                      acceptedStep = importanceSampling();
                   }

                   if(acceptedStep == true){
                       m_acceptedSteps++;
                       if(i >= 0*m_equilibrationFraction*numberOfMetropolisSteps[n]){
                           m_sampler->sample(acceptedStep);

                           //                if ((i%100==0) && (i != 0)){steps ++;}
                           //          Here you should sample the energy (and maybe other things using
                           //          the m_sampler instance of the Sampler class. Make sure, though,
                           //          to only begin sampling after you have let the system equilibrate
                           //          for a while. You may handle this using the fraction of steps which
                           //          are equilibration steps; m_equilibrationFraction.



                           //      }
                           // write only every 10 value

//                           if ((m_numberOfMetropolisSteps > 1e4) && (i%10==0) && (i != 0)){
//                               m_sampler->writeStepToFile(i, steps);
//                               steps++;
//                           }
//                           if ((m_numberOfMetropolisSteps <= 1e4) && (i%1==0) && (i != 0)){
//                               m_sampler->writeStepToFile(i, steps);
//                               steps++;
//                           }
                       }
                   }
                   m_MCstep++;

               }
               m_sampler->computeAverages();
               m_sampler->printOutputToTerminal();
               m_sampler->writeVarToFile();
               cout << "My step index  "<<m_nStepIndex << endl;


               // m_sampler->writeTimeStepToFile();
//               m_sampler->writeAlphaToFile();
//               m_acceptedSteps_ratio = m_acceptedSteps/((double) m_numberOfMetropolisSteps[n]);
               m_acceptedSteps_ratio = m_acceptedSteps/((double)(numberOfMetropolisSteps[n]));
               cout << "Acceptance rate: " << m_acceptedSteps_ratio << endl;
               //        cout << "Acceptance rate importance sampling: " << m_stepImportance/((double) m_numberOfMetropolisSteps) << endl;
               //m_sampler->writeTimeStepToFile();

               m_nStepIndex++;
            }
               m_timestep++;
           }
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


