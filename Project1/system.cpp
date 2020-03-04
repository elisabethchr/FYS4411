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
    random_d = Random::nextDouble();

    s = Random::nextDouble();

    std::vector<double> alphas =  m_waveFunction->getParameters();

    // double oldWaveFunction = m_waveFunction->evaluate(m_particles,alphas[m_sampler->m_stepNumber]); //added m_parameters[m_sampler->m_stepNumber] istead of alpha
    double oldWaveFunction = m_waveFunction->evaluate(m_particles);



    //   double m_position_old, m_position_new;
    //   m_position_old = particles[random_i]->getPosition()[0];
    //   m_position_new = m_position_old + delta*(random_d - 0.5);
    double dx = m_stepLength*(random_d - 0.5);

    m_particles[random_i]->adjustPosition(dx,0);
    double newWaveFunction = m_waveFunction->evaluate(m_particles);

    double ratio = newWaveFunction*newWaveFunction/(oldWaveFunction*oldWaveFunction); //Fix: can simplify expression to save cpu cycles
    //      cout << "Old Wavefunction:" << oldWaveFunction << "    New Wavefunction" << newWaveFunction <<endl;
    if(s<=ratio){
        return true;
    }else{m_particles[random_i]->adjustPosition(-dx,0);
        return false;
    }

    //   Particle->setPosition(m_position_new);

}

void System::runMetropolisSteps(int numberOfMetropolisSteps) {
    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this); //Remove later: (this) points to the system object from which  "runMetropolisSteps" is called.
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    /*
     * InitialPos = GetPosition();
*/

    //    for(const auto i:m_particles){cout << i << endl; }
    for (int i=0; i < numberOfMetropolisSteps; i++) {
        for (int j=0; j < m_numberOfParticles; j++){

            /*Update positions of particles*/
            bool acceptedStep = metropolisStep();
            if(acceptedStep == true){

                m_sampler->sample(acceptedStep);


                /* Here you should sample the energy (and maybe other things using
         * the m_sampler instance of the Sampler class. Make sure, though,
         * to only begin sampling after you have let the system equilibrate
         * for a while. You may handle this using the fraction of steps which
         * are equilibration steps; m_equilibrationFraction.
         */

            }
        }
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






//#include "system.h"
//#include <cassert>
//#include <iostream>
//#include <fstream>
//#include <iomanip>
//#include <armadillo>
//#include "sampler.h"
//#include "particle.h"
//#include "WaveFunctions/wavefunction.h"
//#include "Hamiltonians/hamiltonian.h"
//#include "InitialStates/initialstate.h"
//#include "Math/random.h"

//using namespace std;
//using namespace arma;

//bool System::metropolisStep(int numberOfParticles, int numberOfDimensions, double delta, int i) {
//    /* Perform the actual Metropolis step: Choose a particle at random and
//     * change its position by a random amount, and check if the step is
//     * accepted by the Metropolis test (compare the wave function evaluated
//     * at this new position with the one at the old position).
//     */
//    /*
//     * Look at/draw one particle at a time. Allocate matrices which contain the position of the particle.
//    */
//    double random_d;
//    int random_i;
//    double s;
//    double dx;

//    random_i = Random::nextInt(numberOfParticles);      // Need uniform distribution
//    random_d = Random::nextDouble();
//    s = Random::nextDouble();
//    dx = delta*2*(random_d - 0.5);

//    //    cout << "dx: " << dx << ", delta: " << delta << endl;
//    double oldWaveFunction = m_waveFunction->evaluate(m_particles);
//    double newWaveFunction;
//    double ratio;



//    m_particles[random_i]->adjustPosition(dx,0);
//    newWaveFunction = m_waveFunction->evaluate(m_particles);
//    ratio = newWaveFunction*newWaveFunction/(oldWaveFunction*oldWaveFunction);


//    /*
//    for(int i=0; i<m_particles.size(); i++){
//        oldWaveFunction = m_waveFunction->evaluate(m_particles);

//        random_d = Random::nextDouble();
//        s = Random::nextDouble();
//        dx = delta*2*(random_d - 0.5);

//        m_particles[i]->adjustPosition(dx,0);
//        newWaveFunction = m_waveFunction->evaluate(m_particles);
//        ratio = newWaveFunction*newWaveFunction/(oldWaveFunction*oldWaveFunction);
//    }
//*/
//    /*
//    cout << "oldWaveFunction: " << oldWaveFunction << " ";
//    cout << "newWaveFunction: " << newWaveFunction << "\n";
//    cout << "Ratio = " << ratio << ", s = " << s << "\n";
//*/
//    if(s<ratio){//cout << "True" << endl;
//        return true;
//    }else{m_particles[random_i]->adjustPosition(-dx, 0);
//        return false;
//    }
//}

//void System::runMetropolisSteps(int numberOfMetropolisSteps) {
//    m_particles                 = m_initialState->getParticles();
//    m_sampler                   = new Sampler(this);
//    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
//    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);
//    bool acceptedStep;

//    for (int i=0; i < numberOfMetropolisSteps; i++) {
////        cout << "Step: " << i;
//        /*Update positions of particles*/
//        for(int j=0; j<m_particles.size(); j++){
////            cout << "particle: " << j ;
//            acceptedStep = metropolisStep(m_numberOfParticles, m_numberOfDimensions, m_stepLength, i);
//            // Ensure equilibration of system
//            if(i >= m_equilibrationFraction*numberOfMetropolisSteps){

//                if(acceptedStep == true){
////                    cout << " --> True" << endl;

//                    /* Here you should sample the energy (and maybe other things using
//         * the m_sampler instance of the Sampler class. Make sure, though,
//         * to only begin sampling after you have let the system equilibrate
//         * for a while. You may handle this using the fraction of steps which
//         * are equilibration steps; m_equilibrationFraction.
//         */
//                    m_sampler->sample(acceptedStep);
//                }
//            }
//        }

//    }
//    m_sampler->computeAverages();
//    m_sampler->printOutputToTerminal();
//}


//void System::setNumberOfParticles(int numberOfParticles) {
//    m_numberOfParticles = numberOfParticles;
//}

//void System::setNumberOfDimensions(int numberOfDimensions) {
//    m_numberOfDimensions = numberOfDimensions;
//}

//void System::setStepLength(double stepLength) {
//    assert(stepLength >= 0);
//    m_stepLength = stepLength;
//}

//void System::setEquilibrationFraction(double equilibrationFraction) {
//    assert(equilibrationFraction >= 0);
//    m_equilibrationFraction = equilibrationFraction;
//}

//void System::setHamiltonian(Hamiltonian* hamiltonian) {
//    m_hamiltonian = hamiltonian;
//}

//void System::setWaveFunction(WaveFunction* waveFunction) {
//    m_waveFunction = waveFunction;
//}

//void System::setInitialState(InitialState* initialState) {
//    m_initialState = initialState;
//}


