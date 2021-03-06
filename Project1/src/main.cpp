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


int main() {
    clock_t c_start = clock();

    int n = 10;     // number of values for alpha
    int m = 10;     // number of values for time steps dt
    int k = 10;     // length of array for the number of Metropolis steps
    double alpha_min = 0.2;
    double alpha_max = 0.7;
    double d_alpha = (alpha_max - alpha_min)/n;
    double timestep_max = 1.0;
    double timestep_min = 0.1;
    double dt = (timestep_max - timestep_min)/m;

    double a; // hard-core diamater



    std::vector<double> alpha;      //Variational parameter for the wavefunction
    std::vector<double> beta;       //Variational parameter for elliptical
    std::vector<double> timestep;   //Timestep (delta t) used in Importance sampling
    std::vector<int> MC_cycles;     //Number of Monte Carlo steps in a sample


/////////////////////////////////////////////////////////////////////////////
/// Bools to determine program flow
    bool numeric = false;            //Set type of calculation (numeric or analytic)
    bool bruteForce = true;         //Set type of solver (brute force or importance sampling)

    bool elliptical = true;         //Use elliptical or spherically symmetric potential/wavefunction
    bool Jastrow = true;            //Include or exclude Jastrow factor, i.e. interaction between particles

    bool alphaVec = true;          //Sets alpha to a vector
    bool dtVec = false;             //Sets dt to a vector
    bool nSteps = false;            //sets numberOfSteps to a vector


/////////////////////////////////////////////////////////////////////////////
/// Set parameters

    if (elliptical==false){ beta.push_back(1); }
    else {beta.push_back(sqrt(8)); }

    if (Jastrow==false){a = 0.0; beta.push_back(sqrt(8));}
    else{a=0.0043; beta.push_back(sqrt(8));}

    // set alpha
    if (alphaVec){
        for(int i=0; i<n+1; i++){ alpha.push_back(alpha_min + i*d_alpha);}
    }else{
       alpha.push_back(0.5);          //Set scalar alpha value here
    }

    // set time steps
    if (dtVec){
        for(int i=0; i<m+1; i++){ timestep.push_back(timestep_min + i*dt);}
    }else{
      timestep.push_back(0.01);       //Set scalar dt value here
    }

    // set number of Metropolis steps
    if (nSteps){
        for(int i=0; i<k+1; i++){ MC_cycles.push_back(pow(2, 10+i)); cout << "2^" << 10+i << " = " << MC_cycles[i] << endl;}
    }else{
      MC_cycles.push_back(pow(2,17));  //Set scalar numberOfSteps value here
    }


/////////////////////////////////////////////////////////////////////////////////
/// set main system

    //Set number of bins and maximum radius for one-body density
    int nBins = 500;
    double r_max = 4;
    double r_min = r_max/nBins;


    int numberOfDimensions  = 3;
    int numberOfParticles   = 1;
    double omega            = 1.0;          // Oscillator frequency.
    double gamma            = 2.8483;       // perturbation when elliptical
    double stepLength       = 1.0;          // Metropolis step length.
    double equilibration    = pow(2, 12);          // Amount of the total steps used for equilibration.


    System* system = new System();
    system->setOneBodyDensity           (nBins,r_max, r_min, beta);
    system->setCalculation              (numeric);
    system->setTimeSteps                (timestep);
    system->setSolver                   (bruteForce);
    if(elliptical){
        system->setHamiltonian              (new EllipticalHarmonicOscillator(system, omega));
        system->setWaveFunction             (new EllipticalGaussian(system, alpha, beta, a, Jastrow, gamma));
        cout << "Elliptical"<<endl;
    }else{
        system->setHamiltonian              (new HarmonicOscillator(system, omega));
        system->setWaveFunction             (new SimpleGaussian(system, alpha));
        cout << "Spherical" << endl;
    }
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);
    system->runMetropolisSteps          (MC_cycles);

    clock_t c_end = clock();

    cout << "\n Total CPU-time used: " << (c_end - c_start)/CLOCKS_PER_SEC << " s" << endl;



/////////////////////////////////////////////////////////////////////////////////
/// method for steepest descent

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
