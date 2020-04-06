#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <iomanip>
#include <cstdlib>
#include <algorithm>
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



    // set parameters
    for(int i=0; i<n+1; i++){ alpha.push_back(alpha_min + i*d_alpha); cout << alpha[i] << endl; }
//    alpha.push_back(0.2);

//    beta.push_back(1.0);
    beta.push_back(2.83843);


    // set timesteps
    //    for(int i=0; i<m+1; i++){ tim_omegamestep.push_back(timestep_min + i*dt); cout << timestep[i] << endl;  }
    timestep.push_back(0.4);

    // set number of Monte CarIdet du kjørte programmet forlo cyclesm_omega
    //    for(int i=0; i<k+1; i++){ MC_cycles.push_back(pow(2, 10+i)); cout << "2^" << 10+i << " = " << MC_cycles[i] << endl;}
    MC_cycles.push_back(pow(2, 17));

    // set type of calculation
    bool numeric = true;


    // set type of solver (bruteForce = true -> solver = bruteForce; bruteForce = false -> solver = importance)
    bool bruteForce = false;


    int numberOfDimensions  = 3;
    int numberOfParticles   = 1;
    int numberOfSteps       = (int) 1e4;
    double omega            = 1.0;          // Oscillator frequency.
    double stepLength       = 1.0;          // Metropolis step length.
    double equilibration    = pow(2, 12);          // Amount of the total steps used for equilibration.

    System* system = new System();
    system->setCalculation              (numeric);
    system->setTimeSteps                (timestep);
    system->setSolver                   (bruteForce);
//    system->setHamiltonian              (new EllipticalHarmonicOscillator(system, omega));
//    system->setWaveFunction             (new EllipticalGaussian(system, alpha, beta, a));
    system->setHamiltonian              (new HarmonicOscillator(system, omega));
    system->setWaveFunction             (new SimpleGaussian(system, alpha));
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);
    system->runMetropolisSteps          (MC_cycles);

    clock_t c_end = clock();

    cout << "\n Total CPU-time used: " << (c_end - c_start)/CLOCKS_PER_SEC << " s" << endl;



/*
    // test code, steepest descent
    int nIterations = 10;
    double Energy, EnergyDerivative, alpha_guess, alphagradient, eta, alpha_prev;


    eta = 0.1;
    alpha_guess = 0.1;
    alpha.push_back(alpha_guess);
    System* system = new System();
    system->setCalculation              (numeric);
    system->setTimeSteps                (timestep);
    system->setSolver                   (bruteForce);

    for (int iter=0; iter<nIterations+1; iter++){
        alpha_prev = alpha[0];
        cout << "*************************************************************" << endl;
        cout << "|" << endl;
        cout << "|" << endl;
        cout << "|" << endl;
        cout << "| Iteration no. " << iter << ", alpha = " << alpha_prev << endl;
        cout << "|" << endl;
        cout << "|" << endl;
        cout << "|" << endl;
        cout << "*************************************************************" << endl;

        system->setHamiltonian              (new EllipticalHarmonicOscillator(system, omega));
        system->setWaveFunction             (new EllipticalGaussian(system, alpha, beta, a));
        system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
        system->setEquilibrationFraction    (equilibration);
        system->setStepLength               (stepLength);
        system->runMetropolisSteps          (MC_cycles);

        Energy = system->getEnergy();
        EnergyDerivative = system->getEnergyDerivative();
        alphagradient = EnergyDerivative;
        alpha[0] = alpha_prev - eta*alphagradient;

        ofstream ofile;
        string filename = "data/f/1f_alpha0_";
        string arg1 = to_string(alpha_guess);
        string arg2 = to_string(MC_cycles[0]);
        string arg3 = to_string(eta);
        string arg4 = to_string(stepLength);
        filename.append(arg1);
        filename.append("_learningrate_");
        filename.append(arg3);
        filename.append("_nSteps_");
        filename.append(arg2);
        filename.append("_stepLength_");
        filename.append(arg4);
        filename.append("_.txt");

        if (iter == 0){
            ofile.open(filename, ios::trunc | ios::out);
            ofile << setw(10) << "iteration" << setw(15) << "alpha" <<setw(15) << "<E>" << endl;
        }else{ofile.open(filename, ios::app | ios::out);}

        if (ofile.is_open()){
            ofile << setiosflags(ios::showpoint | ios::uppercase);
            ofile << setw(10) << setprecision(8) << iter;
            ofile << setw(15) << setprecision(8) << alpha_prev;
            ofile << setw(15) << setprecision(8) << Energy << "\n";
            ofile.close();
        }else{
            cout << "Error opening file "<<filename << endl;
        }
    }
*/
    return 0;
}
