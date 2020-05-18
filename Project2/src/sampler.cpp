#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <iomanip>
#include <cstdlib>
#include <algorithm>
#include "sampler.h"
#include "system.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"
#include "Optimizers/optimizer.h"

using namespace std;

Sampler::Sampler(System* system) {
    m_system = system;
    m_stepNumber = 0;
}

void Sampler::setNumberOfMetropolisSteps(int steps) {
    m_numberOfMetropolisSteps = steps;
}


void Sampler::sample(bool acceptedStep) {
    // make sure the sampling variable(s) are initialized at the first step
    if (m_stepNumber == 0) {
        m_nv = m_system->getNumberVisibleNodes();
        m_nh = m_system->getNumberHiddenNodes();
        m_cumulativeEnergy = 0.0;
        m_cumulativeEnergy2 = 0.0;
        m_dPsi = 0.0;
        m_EdPsi = 0.0;
        t_anal = 0;
        m_cumulative_dPsi.zeros(m_nh + m_nv + m_nh*m_nv);
        m_cumulative_EdPsi.zeros(m_nh + m_nv + m_nh*m_nv);
    }

    //calculate cpu-time
    clock_t c_start1 = clock();

    // compute local energy and energy gradient
    arma::vec X = m_system->getWaveFunction()->get_X();
    double localEnergy = m_system->getHamiltonian()->computeLocalEnergy(X);
    arma::vec dPsi = m_system->getHamiltonian()->computeLocalEnergyGradient();

    clock_t c_end1 = clock();
    double t = (c_end1 - c_start1);
    t_num += t;

    m_cumulativeEnergy  += localEnergy;
    m_cumulativeEnergy2 += localEnergy*localEnergy;
    m_cumulative_dPsi   += dPsi;
    m_cumulative_EdPsi  += localEnergy*dPsi;

    m_deltaEnergy = localEnergy;

    m_stepNumber++;

}

void Sampler::computeAverages() {
    // Compute the averages of sampled quantities
    m_energy = m_cumulativeEnergy / m_stepNumber;
    m_energy2 = m_cumulativeEnergy2 / m_stepNumber;
    m_dPsi = m_cumulative_dPsi / m_stepNumber;
    m_EdPsi = m_cumulative_EdPsi / m_stepNumber;

    // compute gradient
    m_gradE = 2*(m_EdPsi - m_energy*m_dPsi);

    // compute variance and error
    m_variance = (m_energy2 - m_energy*m_energy) / m_stepNumber;
    m_error = pow(abs(m_variance), 0.5);

}

/* Optimize weights */
void Sampler::optimizeWeights(){
    m_system->getOptimizer()->computeWeights(m_gradE);
}

/* Display variables to terminal */
void Sampler::printOutputToTerminal() {
    int     nv = m_system->getNumberVisibleNodes();
    int     nh = m_system->getNumberHiddenNodes();
    int     ms = m_system->getNumberOfMetropolisSteps();
    int     bm = m_system->getNumberRBMcycles();
    double  ef = m_system->getEquilibrationFraction();

    cout << endl;
    cout << "  -- System info -- " << endl;
    cout << " Number of visible nodes  : " << nv << endl;
    cout << " Number of hidden nodes : " << nh << endl;
    cout << " Number of RBM cycles run : " << bm << endl;
    cout << " Number of Metropolis steps run : 10^" << std::log10(ms) << endl;
    cout << " Number of equilibration steps  : 10^" << std::log10(std::round(ef)) << endl;
    cout << "Learning rate: " << m_system->getLearningRate() << endl;
    cout << endl;
    cout << "  -- Results -- " << endl;
    cout << " Energy : " << m_energy << endl;
    cout << "Variance: " << m_variance << endl;
    cout << "Error: " << m_error << endl;
}

/* Write values to file for each optimization cycle */
void Sampler::writeToFile(int step, int steps, string filename){

    int nv = m_system->getNumberVisibleNodes();
    int nh = m_system->getNumberHiddenNodes();
    int nP = m_system->getNumberParticles();
    int nD = m_system->getNumberDimensions();
    int h  = m_system->getStepLength();
    double dt = m_system->getTimeSteps()[0];
    int nSteps = m_system->getNumberOfMetropolisSteps();
    int nRBMcycles = m_system->getNumberRBMcycles();
    double eta    = m_system->getLearningRate();
    bool gauss = m_system->getWaveFunction()->getGaussian();
    bool solv = m_system->getSolver();

    string solver;
    if(solv==true){ solver = "bruteForce"; }
    else if (solv==false){ solver = "importance"; }

    string gaussian;
    if(gauss == true){ gaussian = "normal"; }
    else if(gauss == false){ gaussian = "uniform"; }

    ofstream ofile;
    string arg1 = to_string(int(nv));
    string arg2 = to_string(int(nh));
    string arg3 = to_string(int(nSteps));
    string arg4 = to_string(int(nP));
    string arg5 = to_string(int(nD));
    string arg6 = to_string(float(h));
    string arg7 = to_string(float(eta));
    string arg8 = to_string(int(nRBMcycles));
    string arg9 = to_string(double(dt));
    filename.append(solver);
    filename.append("_");
    filename.append(gaussian);
    filename.append("_nv_");
    filename.append(arg1);
    filename.append("_nh_");
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
    filename.append(arg7);
    filename.append("_RBMcycles_");
    filename.append(arg8);
    if(solv==false){
        filename.append("_dt_");
        filename.append(arg9);
    }
    filename.append("_.txt");

    if (steps == 0){
        ofile.open(filename, ios::trunc | ios::out);
        ofile << setw(10) << "steps" <<setw(15) << "Energy" << setw(15) << "Variance" << setw(15) << "Error" << endl;
    }
    else{ofile.open(filename, ios::app | ios::out);}

    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(10) << setprecision(8) << step;
    ofile << setw(15) << setprecision(8) << m_energy;
    ofile << setw(15) << setprecision(8) << m_variance;
    ofile << setw(15) << setprecision(8) << m_error << endl;
    ofile.close();
}


/* Write sampled quantities for each Metropolis step to file */
void Sampler::writeStepToFile(int step, int steps, string filename){
    int nv = m_system->getNumberVisibleNodes();
    int nh = m_system->getNumberHiddenNodes();
    int nP = m_system->getNumberParticles();
    int nD = m_system->getNumberDimensions();
    int h  = m_system->getStepLength();
    double dt = m_system->getTimeSteps()[0];
    double eta    = m_system->getLearningRate();
    double nSteps = m_system->getNumberOfMetropolisSteps();
    double nRBMcycles = m_system->getNumberRBMcycles();
    bool gauss = m_system->getWaveFunction()->getGaussian();
    bool solv = m_system->getSolver();

    string solver;
    if(solv==true){ solver = "bruteForce"; }
    else if (solv==false){ solver = "importance"; }

    string gaussian;
    if(gauss == true){ gaussian = "normal"; }
    else if(gauss == false){ gaussian = "uniform"; }

    ofstream ofile;
    string arg1 = to_string(int(nv));
    string arg2 = to_string(int(nh));
    string arg3 = to_string(int(nSteps));
    string arg4 = to_string(int(nP));
    string arg5 = to_string(int(nD));
    string arg6 = to_string(float(h));
    string arg7 = to_string(float(eta));
    string arg8 = to_string(int(nRBMcycles));
    string arg9 = to_string(double(dt));
    filename.append(solver);
    filename.append("_");
    filename.append(gaussian);
    filename.append("_nv_");
    filename.append(arg1);
    filename.append("_nh_");
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
    filename.append(arg7);
    filename.append("_RBMcycles_");
    filename.append(arg8);
    if(solv==false){
        filename.append("_dt_");
        filename.append(arg9);
    }
    filename.append("_.txt");

    if (steps == 0){
        ofile.open(filename, ios::trunc | ios::out);
        ofile << setw(10) << "steps" <<setw(15) << "Energy" << endl;
        cout << "step = " << step << ", m_deltaEnergy = " << m_deltaEnergy << endl;
    }
    else{ofile.open(filename, ios::app | ios::out);}

    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(10) << setprecision(8) << step;
    ofile << setw(15) << setprecision(8) << m_deltaEnergy << endl;
    ofile.close();
}



/* Write averaged quantities for each value of time step (dt), when running importance sampling, to file */
void Sampler::writeTimeStepToFile(string filename){

    int i = m_system->getTimeStepIndex();
    int step = m_system->getSampleStep();
    double nParticles = m_system->getNumberParticles();
    double nDim = m_system->getNumberDimensions();
    double nSteps = m_system->getNumberOfMetropolisSteps();
    double timestep = m_system->getTimeSteps()[i];
    double stepLength = m_system->getStepLength();

    ofstream ofile;
    string arg1 = to_string(int(nParticles));
    string arg2 = to_string(int(nDim));
    string arg3 = to_string(int(nSteps));
    string arg4 = to_string(stepLength);
    filename.append(arg1);
    filename.append("_nDim_");
    filename.append(arg2);
    filename.append("_nSteps_");
    filename.append(arg3);
    filename.append("_stepLength_");
    filename.append(arg4);
    filename.append("_");
    filename.append(".txt");
    if (i == 0){
        ofile.open(filename, ios::trunc | ios::out);
        ofile << setw(7) << "step" << setw(15) << "Energy" << setw(15) << "Variance" << setw(15) << "Error" << setw(15) << "ratio_steps" << "\n";
    }
    else{ofile.open(filename, ios::app | ios::out);}

    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(10) << setprecision(8) << timestep;
    ofile << setw(15) << setprecision(8) << m_cumulativeEnergy / ((double) m_stepNumber);
    ofile << setw(15) << setprecision(8) << m_variance;
    ofile << setw(15) << setprecision(8) << m_error;
    ofile << setw(15) << setprecision(8) << m_stepNumber / ((double) m_system->getNumberOfMetropolisSteps()) << "\n";
    ofile.close();

}
