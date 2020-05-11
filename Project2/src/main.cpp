#include <iostream>
#include <ctime>
#include <armadillo>
#include "system.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/neuralquantumstate.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/initalizestate.h"
#include "Optimizers/optimizer.h"
#include "Optimizers/stochasticgradientdescent.h"
#include "Math/random.h"
#include "sampler.h"

using std::cout;
using std::endl;
using namespace std;
using namespace arma;


int main() {
    clock_t c_start = clock();

    int m = 10;     // number of values for time steps dt
    int k = 10;     // length of array for the number of Metropolis steps
    int MC = pow(2, 17);
    double eq = round(log2(0.1*MC));
    double timestep_max = 1.0;
    double timestep_min = 0.1;
    double dt = (timestep_max - timestep_min)/m;

    std::vector<int> MC_cycles;     //Number of Monte Carlo steps in a sample
    std::vector<double> timestep;   //Timestep (delta t) used in Importance sampling


/////////////////////////////////////////////////////////////////////////////
/// Bools to determine program flow
    bool bruteForce = true;         // set type of solver (brute force or importance sampling)
    bool interaction = false;       // if interaction should be included or not

    bool alphaVec = false;          // sets alpha to a vector
    bool dtVec = false;             // sets dt to a vector
    bool nSteps = false;            // sets numberOfSteps to a vector


/////////////////////////////////////////////////////////////////////////////
/// Set parameters

    // set time steps
    if (dtVec){
        for(int i=0; i<m+1; i++){ timestep.push_back(timestep_min + i*dt);}
    }else{
      timestep.push_back(0.45);       //Set scalar dt value here
    }

    // set number of Metropolis steps
    if (nSteps){
        for(int i=0; i<k+1; i++){ MC_cycles.push_back(pow(2, 10+i)); cout << "2^" << 10+i << " = " << MC_cycles[i] << endl;}
    }else{
      MC_cycles.push_back(pow(2,17));  //Set scalar numberOfSteps value here
    }


/////////////////////////////////////////////////////////////////////////////////
/// set main system
    int P = 1;      // number of particles
    int D = 1;      // number of dimensions

    int RBM_cycles = 20;       // number of sampling cycles
    int nVisible  = P*D;    // number of visible nodes
    int nHidden   = 2;      // number of hidden nodes
    double omega            = 1.0;          // Oscillator frequency
    double stepLength       = 1.0;          // Metropolis step length
    double sigma            = 1.0;          // error of Gaussian distribution
    double equilibration    = pow(2, eq);   // Amount of the total steps used for equilibration.
    cout << "Equilibration: 2^" << eq << " = " << equilibration << endl;
    double eta            = 0.1;          // Learning rate

    System* system = new System();
//    system->setTimeSteps                (timestep);
    system->setSolver                   (bruteForce);
    system->setHamiltonian              (new HarmonicOscillator(system, omega, interaction));
    system->setWaveFunction             (new NeuralQuantumState(system, nHidden, nVisible, P, D, sigma));
    system->setOptimizer                (new StochasticGradientDescent(system, eta));
    //system->setSampler                  (new Sampler(system));
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);
    system->runMetropolisSteps          (RBM_cycles, MC_cycles);

    clock_t c_end = clock();

    cout << "\n Total CPU-time used: " << (c_end - c_start)/CLOCKS_PER_SEC << " s" << endl;

    return 0;
}
