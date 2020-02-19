#include <iostream>
#include <ctime>
#include <armadillo>
#include "system.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "Math/random.h"

using namespace std;
using namespace arma;


/*
 * See hamiltonian.h for explanation on class <name1> : public <name2> (inheritages of member calls and member variables)
*/

int main() {
    clock_t c_start = clock();

    int numberOfDimensions  = 1;
    int numberOfParticles   = 1;
    int numberOfSteps       = (int) 2;
    double omega            = 1.0;          // Oscillator frequency.
    double alpha            = 0.5;          // Variational parameter.
    double stepLength       = 0.1;          // Metropolis step length.
    double equilibration    = 0.1;          // Amount of the total steps used
    // for equilibration.

    System* system = new System();
    system->setHamiltonian              (new HarmonicOscillator(system, omega));
    system->setWaveFunction             (new SimpleGaussian(system, alpha));
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);
    system->runMetropolisSteps          (numberOfSteps);

    clock_t c_end = clock();

    cout << "CPU-time used: " << (c_end - c_start)/1000. << "ms" << endl;

    return 0;
}
