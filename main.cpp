#include <iostream>
#include "system.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "Math/random.h"

#include <fstream>

using namespace std;


//int n = 4;                  //Number of elements in nParticles
//int nParticles[4] = {1,10,20,100};

int main() {



//ofstream data;
//data.open("data4days.txt");
//data << "#CPU run time for varying particle numbers, 100 steps"<<endl;
//data << "n t[ms]" <<endl;
////data << "Number of particles    CPU time"<<endl;


bool importance = false; //Set to true for Importance, false for brute force Metropolis

int alphalength = 40;
double alpha_max = 1;
double alpha_min = 0.2;
double dalpha = (alpha_max-0.2)/alphalength;
//vec<double> alphas =

//for(int i = 0;i<n;i++){
//       clock_t c_start = clock();
   for(int j = 0; j<alphalength; j++){

        int numberOfDimensions  = 1;
//        int numberOfParticles   = nParticles[i];
        int numberOfParticles   = 1;
        int numberOfSteps       = (int) 1e7;
        double omega            = 1.0;          // Oscillator frequency.
        double alpha            = alpha_min+j*dalpha;          // Variational parameter.
        double stepLength       = 0.1;          // Metropolis step length.
        double equilibration    = 0.1;          // Amount of the total steps used
        // for equilibration.

        System* system = new System();
        system->setHamiltonian              (new HarmonicOscillator(system, omega));
        system->setWaveFunction             (new SimpleGaussian(system, alpha));
        system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
        system->setEquilibrationFraction    (equilibration);
        system->setStepLength               (stepLength);
        system->runMetropolisSteps          (numberOfSteps,importance);
//        cout << "Analytical Energy : "<< numberOfDimensions/2*numberOfParticles*omega <<endl;
//        clock_t c_end = clock();

//        cout << "CPU-time used: " << (c_end - c_start)/1000. << "ms" << endl;
//        data <<nParticles[i]<<"         " << (c_end-c_start)/1000. <<endl;
    }
//}
//data.close();

    return 0;
}
