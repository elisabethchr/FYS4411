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
    //    clock_t c_start = clock();

    // setting variational parameter alpha
    int n = 30;
    double alpha_min = 0.2;
    double alpha_max = 0.7;
    double d_alpha = (alpha_max - alpha_min)/n;
    std::vector<double> alpha;

    int m = 100;
    double timestep_max = 0.1;
    double timestep_min = 0.0001;
    double dt = (timestep_max - timestep_min)/m;
    std::vector<double> timestep;

    int l = 10;
    int nSteps_max = 1e5;
    int nSteps_min = 1e2;
    int dn = (nSteps_max - nSteps_min)/l;
    std::vector<int> numberOfSteps;


    bool numeric = true;
    bool bruteForce = false;

    bool alphaVec = false;       //Sets alpha to a vector
    bool dtVec = false;         //Sets dt to a vector
    bool nSteps = true;         //sets numberOfSteps to a vector

    if (alphaVec){
        for(int i=0; i<n+1; i++){ alpha.push_back(alpha_min + i*d_alpha);}
    }else{
       alpha.push_back(0.3);          //Set scalar alpha value here
    }
    if (dtVec){
        for(int i=0; i<m+1; i++){ timestep.push_back(timestep_min + i*dt);}
    }else{
      timestep.push_back(0.1);       //Set scalar dt value here
    }
    if (nSteps){
        for(int i=0; i<l+1; i++){ numberOfSteps.push_back(nSteps_min + i*dn); cout << numberOfSteps[i]<<endl;}
    }else{
      numberOfSteps.push_back(1e4);  //Set scalar numberOfSteps value here
    }


    int numberOfDimensions  = 1;
    int numberOfParticles   = 1;
 //   int numberOfSteps       = (int) 1e4;
    double omega            = 1.0;          // Oscillator frequency.
    double stepLength       = 0.5;          // Metropolis step length.
    double equilibration    = 0.05;          // Amount of the total steps used for equilibration.
    System* system = new System();
    system->setCalculation              (numeric);
    system->setTimeSteps                (timestep);
    system->setHamiltonian              (new HarmonicOscillator(system, omega));
    system->setWaveFunction             (new SimpleGaussian(system, alpha));
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);
    system->runMetropolisSteps          (numberOfSteps, bruteForce);



//    clock_t c_end = clock();

//    cout << "Total CPU-time used: " << (c_end - c_start)/1000. << "ms" << endl;

    return 0;
}
