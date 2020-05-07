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
#include "Math/random.h"

using std::cout;
using std::endl;
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

    std::vector<double> alpha;      //Variational parameter for the wavefunction
    std::vector<double> beta;       //Variational parameter for elliptical
    std::vector<double> timestep;   //Timestep (delta t) used in Importance sampling
    std::vector<int> MC_cycles;     //Number of Monte Carlo steps in a sample


/////////////////////////////////////////////////////////////////////////////
/// Bools to determine program flow
    bool bruteForce = false;         //Set type of solver (brute force or importance sampling)

    bool alphaVec = false;          //Sets alpha to a vector
    bool dtVec = false;             //Sets dt to a vector
    bool nSteps = false;            //sets numberOfSteps to a vector


/////////////////////////////////////////////////////////////////////////////
/// Set parameters
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
      timestep.push_back(0.45);       //Set scalar dt value here
    }

    // set number of Metropolis steps
    if (nSteps){
        for(int i=0; i<k+1; i++){ MC_cycles.push_back(pow(2, 10+i)); cout << "2^" << 10+i << " = " << MC_cycles[i] << endl;}
    }else{
      MC_cycles.push_back(pow(2,20));  //Set scalar numberOfSteps value here
    }


/////////////////////////////////////////////////////////////////////////////////
/// set main system
    int P = 2;      // number of particles
    int D = 2;      // number of dimensions

    int nVisible  = P*D;    // number of visible nodes
    int nHidden   = 2;     // number of hidden nodes
    double omega            = 1.0;          // Oscillator frequency.
    double stepLength       = 1.0;          // Metropolis step length.
    double equilibration    = pow(2, 12);          // Amount of the total steps used for equilibration.


    System* system = new System();
    system->setTimeSteps                (timestep);
    system->setSolver                   (bruteForce);
    system->setHamiltonian              (new HarmonicOscillator(system, omega));
    system->setWaveFunction             (new NeuralQuantumState(system));
    system->setInitialState             (new InitializeState(system, nHidden, nVisible, D));
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);
    system->runMetropolisSteps          (MC_cycles);

    clock_t c_end = clock();

    cout << "\n Total CPU-time used: " << (c_end - c_start)/CLOCKS_PER_SEC << " s" << endl;

    return 0;
}