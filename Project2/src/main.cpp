#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <armadillo>
#include "system.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/neuralquantumstate.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
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
    int MC = 19;    // log2 of number of Monte Carlo cycles
    int MC_pow2 = pow(2, MC);   // total number of Monte Carlo cycles
    double eq = round(log2(0.1*MC_pow2));
    double timestep_max = 1.0;
    double timestep_min = 0.1;
    double dt = (timestep_max - timestep_min)/m;

    std::vector<int> MC_cycles;     //Number of Monte Carlo steps in a sample


    /////////////////////////////////////////////////////////////////////////////
    /// Options to determine program flow

    // Only one solver can be true at a time (!)
    bool bruteForce = false;             // set type of solver to brute force
    bool importance = false;             // set type of solver to importance sampling
    bool gibbs = true;                   // set type of solver to Gibbs sampling

    bool interaction = true;           // if interaction should be included or not
    bool GaussianDistribution = true;   // distribution of weights (either uniform or gaussian)

    bool dtVec = false;             // sets dt to a vector
    bool nSteps = false;            // sets numberOfSteps to a vector


    /////////////////////////////////////////////////////////////////////////////
    /// Set parameters

    // set number of Metropolis steps
    if (nSteps){
        for(int i=0; i<k+1; i++){ MC_cycles.push_back(pow(2, 10+i)); cout << "2^" << 10+i << " = " << MC_cycles[i] << endl;}
    }else{
        MC_cycles.push_back(pow(2,MC));  //Set scalar numberOfSteps value here
    }

    /////////////////////////////////////////////////////////////////////////////////
    /// set main system
    int P = 2;      // number of particles
    int D = 1;      // number of dimensions

    int RBM_cycles = 100;    // number of sampling cycles
    int nVisible   = P*D;    // number of visible nodes
    int nHidden    = 2;      // number of hidden nodes
    double initialization = 0.01;    // initialization interval for weights (if uniform distribution -> uniform(-init, +init), if normal distribution -> normal(sigma = init))
    double omega            = 1.0;          // Oscillator frequency
    double stepLength       = 2.0;          // Metropolis step length
    double sigma            = 1.0;          // error of Gaussian distribution
    double equilibration    = pow(2, eq);   // Amount of the total steps used for equilibration.
    double eta              = 1.0;          // Learning rate
    arma::vec sigmas = {1.0};//{1.0, 0.75, 0.5, sqrt(0.5), 0.25};//{1.0};
    arma::vec rates = {1.0};//{1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 1.1, 1.2};
    arma::vec nh = {4};//{2, 4, 6, 8, 10, 12, 14, 16, 18, 20};
    arma::vec timestep = {1.35};    // timestep (dt) used in importance sampling

    /////////////////////////////////////////////////////////////////////////////////
    /// the rest is automatic (unless any changes to filenames are wanted)

    string filename;
    string arg0 = to_string(initialization);
    if (nh.size() > 1){
        filename = "../data/g/hiddenNodes/2g-Initialization-";
        filename.append(arg0);
        filename.append("_hiddenNodes_");
    }
    else if(rates.size() > 1){
        filename = "../data/g/learningRate/2g-Initialization-";
        filename.append(arg0);
        filename.append("_learningRate_");
    }
    else if(sigmas.size() > 1){
        filename = "../data/g/variance/2g-Initialization-";
        filename.append(arg0);
        filename.append("_sigma_");
    }

    string solver;
    if ((bruteForce==true && importance == true) || (bruteForce==true && gibbs==true) || (importance==true && gibbs==true)){
        cout << "Only one solver can be true at a time!" << endl;
        exit (EXIT_FAILURE);
    }
    else if (bruteForce==false && importance==false && gibbs==false){
        cout << "You must set a solver!" << endl;
        exit (EXIT_FAILURE);
    }

    if (bruteForce==true){ solver = "bruteForce"; }
    else if (importance==true){ solver = "importance"; }
    else if (gibbs == true){ solver = "gibbs"; }

    string gaussian;
    if(GaussianDistribution == true){ gaussian = "normal"; }
    else if(GaussianDistribution == false){ gaussian = "uniform"; }

    // write to file
//    ofstream ofile;
    string arg1 = to_string(nVisible);
    string arg2;        // number hidden nodes (can be variable)
    string arg3 = to_string(MC_cycles[0]);
    string arg4 = to_string(P);
    string arg5 = to_string(D);
    string arg6 = to_string(stepLength);
    string arg7;        // learning rate (can be variable)
    string arg8 = to_string(RBM_cycles);
    string arg9;
    //    string arg4 = to_string(double(dt));
    filename.append(solver);
    filename.append("_");
    filename.append(gaussian);
    filename.append("_nv_");
    filename.append(arg1);
    filename.append("_nh_");
    if (nh.size() > 1){ arg2 = "x"; }
    else if (nh.size() == 1){ arg2 = to_string(int(nh[0])); }
    filename.append(arg2);
    filename.append("_nMCsteps_");
    filename.append(arg3);
    filename.append("_nPart_");
    filename.append(arg4);
    filename.append("_nDim_");
    filename.append(arg5);
    filename.append("_stepL_");
    filename.append(arg6);
    filename.append("_eta_");
    if (rates.size() > 1){ arg7 = "x"; }
    else if (rates.size() == 1){ arg7 = to_string(rates[0]); }
    filename.append(arg7);
    filename.append("_RBMcycles_");
    filename.append(arg8);
    if (solver=="gibbs"){
        filename.append("_sigma_");
        if (sigmas.size() > 1){arg9 = "x"; }
        else if (sigmas.size() == 1){ arg9 = to_string(sigmas[0]); }
    }
    filename.append("_.txt");

    // define length of loop
    int N;
    if (nh.size() > 1){ N = nh.size(); }
    else if (rates.size() > 1){ N = rates.size(); }
    else if (sigmas.size() > 1){ N = sigmas.size(); }
    int i = 0;
//    for (int i=0; i<N; i++){
        // set variables for learning rate and number hidden nodes
        int a=0;   // number hidden nodes
        double b=0.0;   // learning rate
        double c=0.0;   // variance (gaussian error)
        if (nh.size() > 1){ a = nh[i]; }
        else{ a = nh[0]; }
        if (rates.size() > 1){ b = rates[i]; }
        else { b = rates[0]; }
        if (sigmas.size() > 1){ c = sigmas[i]; }
        else { c = sigmas[0]; }

        System* system = new System();
        system->setTimeSteps                (timestep);
        system->setSolver                   (solver);
        system->setHamiltonian              (new HarmonicOscillator(system, omega, interaction));
        system->setWaveFunction             (new NeuralQuantumState(system, a, nVisible, P, D, c, GaussianDistribution, initialization));
        system->setOptimizer                (new StochasticGradientDescent(system, b));
        system->setEquilibrationFraction    (equilibration);
        system->setStepLength               (stepLength);
        system->runMetropolisSteps          (RBM_cycles, MC_cycles);

        double m_energy = system->getSampler()->getEnergy();
        double m_error = system->getSampler()->getError();

//        if (i == 0){
//            ofile.open(filename, ios::trunc | ios::out);
//            ofile << setw(10) << "steps" <<setw(15) << "Energy" << setw(15) << "Error" << endl;
//        }
//        else{ofile.open(filename, ios::app | ios::out);}

//        if (nh.size()>1){ cout << "nh = " << a << ", E_L = " << m_energy << endl; }
//        else if (rates.size()>1){ cout << "eta = " << b << ", E_L = " << m_energy << endl; }
//        else if (sigmas.size()>1){ cout << "sigma = " << c << ", E_L = " << m_energy << endl; }

//        ofile << setiosflags(ios::showpoint | ios::uppercase);
//        if (nh.size()>1){ ofile << setw(10) << setprecision(8) << nh[i]; }
//        else if (rates.size()>1){ ofile << setw(10) << setprecision(10) << rates[i]; }
//        else if (sigmas.size()>1){ ofile << setw(10) << setprecision(10) << sigmas[i]; }
//        ofile << setw(15) << setprecision(8) << m_energy;
//        ofile << setw(15) << setprecision(8) << m_error << endl;
//        ofile.close();
//    }
    clock_t c_end = clock();

    cout << "\n Total CPU-time used: " << (c_end - c_start)/CLOCKS_PER_SEC << " s" << endl;

    return 0;
}
