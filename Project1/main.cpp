#include <iostream>
#include <ctime>
#include <armadillo>
#include "system.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "WaveFunctions/ellipticalgaussian.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "Hamiltonians/elliptical_harmonicoscillator.h"
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

    // setting variational parameter alpha
    int n = 10;
    int m = 10;
    int k = 10;
    double alpha_min = 0.2;
    double alpha_max = 0.7;
    double d_alpha = (alpha_max - alpha_min)/n;
    double timestep_max = 1.0;
    double timestep_min = 0.1;
    double dt = (timestep_max - timestep_min)/m;
    double a = 0.0043;
//    double a = 0.0;

    std::vector<double> alpha;
    std::vector<double> beta;
    std::vector<double> timestep;
    std::vector<int> MC_cycles;
    bool Jastrow = true;



    // set parameters
    for(int i=0; i<n+1; i++){ alpha.push_back(alpha_min + i*d_alpha); cout << alpha[i] << endl; }
//    alpha.push_back(0.5);

    beta.push_back(2.83843);
//    beta.push_back(1);


    // set timesteps
//    for(int i=0; i<m+1; i++){ tim_omegamestep.push_back(timestep_min + i*dt); cout << timestep[i] << endl;  }
    timestep.push_back(0.01);


    // set number of Monte Carlo cyclesm_omega
//    for(int i=0; i<k+1; i++){ MC_cycles.push_back(pow(2, 10+i)); cout << "2^" << 10+i << " = " << MC_cycles[i] << endl;}
    MC_cycles.push_back(pow(2, 15));
//     MC_cycles.push_back(1e4);


    // set type of calculation
    bool numeric = false;

    // set type of solver (bruteForce = true -> solver = bruteForce; bruteForce = false -> solver = importance)
    bool bruteForce = true;

    //Set number of bins and maximum radius for one-body density
    int nBins = 100;
    double r_max = 3;


    int numberOfDimensions  = 3;
    int numberOfParticles   = 10;
    int numberOfSteps       = (int) 1e6;
    double omega            = 1.0;          // Oscillator frequency.
    double stepLength       = 0.1;          // Metropolis step length.
    double equilibration    = pow(2, 12);          // Amount of the total steps used for equilibration.



//        }
    System* system = new System();
    system->setOneBodyDensity           (nBins,r_max);
    system->setCalculation              (numeric);
    system->setTimeSteps                (timestep);
    system->setSolver                   (bruteForce);
    system->setHamiltonian              (new EllipticalHarmonicOscillator(system, omega));
    system->setWaveFunction             (new EllipticalGaussian(system, alpha, beta, a, Jastrow));
//    system->setHamiltonian              (new HarmonicOscillator(system, omega));
//    system->setWaveFunction             (new SimpleGaussian(system, alpha));
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);
    system->runMetropolisSteps          (MC_cycles);

    clock_t c_end = clock();

    cout << "\n Total CPU-time used: " << (c_end - c_start)/CLOCKS_PER_SEC << " s" << endl;




// test code, steepest descent
//int nIterations = 50;
//double Energy, EnergyDerivative, alpha_guess, alphagradient, eta;


//eta = 0.1;
//alpha_guess = 0.2;
//alpha.push_back(alpha_guess);
//System* system = new System();
//system->setCalculation              (numeric);
//system->setTimeSteps                (timestep);
//system->setSolver                   (bruteForce);

//    for (int iter=0; iter<nIterations; iter++){
//        cout << "*************************************************************" << endl;
//        cout << "|" << endl;
//        cout << "|" << endl;
//        cout << "|" << endl;
//        cout << "| Iteration no. " << iter << ", alpha = " << alpha[0] << endl;
//        cout << "|" << endl;
//        cout << "|" << endl;
//        cout << "|" << endl;
//        cout << "*************************************************************" << endl;

//        system->setHamiltonian              (new EllipticalHarmonicOscillator(system, omega));
//        system->setWaveFunction             (new EllipticalGaussian(system, alpha, beta, a));
//    //    system->setHamiltonian              (new HarmonicOscillator(system, omega));
//    //    system->setWaveFunction             (new SimpleGaussian(system, alpha));
//        system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
//        system->setEquilibrationFraction    (equilibration);
//        system->setStepLength               (stepLength);
//        system->runMetropolisSteps          (MC_cycles);

//        Energy = system->getEnergy();
//        EnergyDerivative = system->getEnergyDerivative();
//        alphagradient = EnergyDerivative;
//        alpha[0] -= eta*alphagradient;


//    }




    return 0;
}
