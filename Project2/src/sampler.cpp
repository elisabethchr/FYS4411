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
#include "particle.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"

using namespace std;

Sampler::Sampler(System* system) {
    m_system = system;
    m_stepNumber = 0;

}

void Sampler::setNumberOfMetropolisSteps(int steps) {
    m_numberOfMetropolisSteps = steps;
}

void Sampler::setOneBodyDensity(int nBins, double r_max, double r_min){
    m_nBins = nBins;
    m_Rmax = r_max;
    m_Rmin = r_min;
    double bin_size = (r_max-r_min)/(double)nBins;
        for(int j = 0; j<nBins; j++){
            m_radii.push_back(r_min+j*bin_size);
            m_OneBodyBin.push_back(0);
        }

}

void Sampler::sample(bool acceptedStep) {
    // Make sure the sampling variable(s) are initialized at the first step.
    if (m_stepNumber == 0) {
        m_cumulativeEnergy = 0;
        m_cumulativeEnergyAnalytic = 0;
        t_num = 0;
        t_anal = 0;
    }

    bool calc = m_system->getCalculation();

    //calculate cpu-time
    clock_t c_start1 = clock();

    double localEnergy = m_system->getHamiltonian()->
            computeLocalEnergy(m_system->getParticles(), calc);
    m_system->setEnergy(localEnergy);


    clock_t c_end1 = clock();
    double t = (c_end1 - c_start1);
    t_num += t;


    m_cumulativeEnergy  += localEnergy;
    m_cumulativeEnergy2 += localEnergy*localEnergy;

    m_derivativePsi_alpha = m_system->getWaveFunction()->computeDerivativePsi_alpha(m_system->getParticles());
    m_deltaEnergy = localEnergy;
    m_cumulativeDeltaPsi += m_derivativePsi_alpha;
    m_cumulativeDerivativePsiE += m_derivativePsi_alpha*m_deltaEnergy;

//    cout << "m_cumulativeEnergy_sample = " << m_cumulativeEnergy << endl;
//    cout << "m_stepNumber_sample = " << m_stepNumber << endl;
//    cout << "m_energy_sample = " << m_cumulativeEnergy/((double)m_stepNumber) << endl;

    m_stepNumber++;
}

void Sampler::sampleOneBodyDensity(){

    //////////////////////////////////////////////////////////////////////
    /// Sample particle positions for one-body density
    ///
     int n_p =m_system->getNumberOfParticles();
     int dim =m_system->getNumberOfDimensions();

     double r_part;
     double dvolume;
     double bin_size = m_Rmax/m_nBins;
     double beta = m_system->getBeta()[0];
     double pi = 3.14159;

//     m_oneBodyBin(n_p,n_bins);

     for(int i = 0; i<n_p; i++){
         std::vector<class Particle*> particles = m_system->getParticles();
             r_part = 0;
             for(int j=0; j<dim; j++){
                 std::vector<double> position = particles[i]->getPosition();
                 if(j<2){
                  r_part +=position[j]*position[j];
                 }else{
                  r_part += beta*position[j]*position[j];
                 }

             }
             r_part = sqrt(r_part);
             if(r_part<=m_Rmax){
                 int bin_i = (int)floor(r_part/bin_size);
//                 dvolume = (4*pi*(pow(m_radii[bin_i],2))*bin_size)/sqrt(beta);
                  dvolume = (4*pi*(pow(m_radii[bin_i],3)-pow(m_radii[bin_i]-bin_size,3)))/(3*sqrt(beta));

                 m_OneBodyBin[bin_i]+= 1/(n_p*dvolume);
             }
     }

}

void Sampler::printOutputToTerminal() {
    int     np = m_system->getNumberOfParticles();
    int     nd = m_system->getNumberOfDimensions();
    int     ms = m_system->getNumberOfMetropolisSteps();
    int     p  = m_system->getWaveFunction()->getNumberOfParameters();
    double  ef = m_system->getEquilibrationFraction();
    std::vector<double> pa = m_system->getWaveFunction()->getParameters();

    cout << endl;
    cout << "  -- System info -- " << endl;
    cout << " Number of particles  : " << np << endl;
    cout << " Number of dimensions : " << nd << endl;
    cout << " Number of Metropolis steps run : 10^" << std::log10(ms) << endl;
    cout << " Number of equilibration steps  : 10^" << std::log10(std::round(ms*ef)) << endl;
    cout << endl;
    cout << "  -- Wave function parameters -- " << endl;
    cout << " Number of parameters : " << p << endl;
    cout << endl;
    cout << "  -- Results -- " << endl;
    cout << " Energy : " << m_energy << endl;
    cout << "Variance: " << m_variance << endl;
    cout << "Error: " << m_error << endl;
}

void Sampler::computeAverages() {
    // Compute the averages of the sampled quantities

    int m_metropolisStep = m_system->getMetropolisStep();

    // test to make sure m_energy != nan
    if(m_stepNumber==0){m_energy = m_cumulativeEnergy; cout << "m_stepNumber = 0" << endl; }
    else {m_energy = m_cumulativeEnergy / ((double) m_stepNumber);}

    m_energy2 = m_cumulativeEnergy2 / ((double )m_stepNumber);
    m_derivativePsiE = m_cumulativeDerivativePsiE / ((double) m_stepNumber);
    m_deltaPsi = m_cumulativeDeltaPsi / ((double) m_stepNumber);
    m_derivativeE = 2*(m_derivativePsiE - m_deltaPsi*m_energy);

    m_variance = (m_energy2 - m_energy*m_energy) / ((double) m_stepNumber);
    m_error = pow(abs(m_variance), 0.5);
}


void Sampler::writeTotalToFile(){
    int nParticles = m_system->getNumberOfParticles();

    ofstream ofile;
    string filename = "10steps.txt";

    if (nParticles == 1){
        ofile.open(filename, ios::trunc | ios::out);
        ofile << setw(10) << "N_{particles}" <<setw(15) << "t_{num}" << setw(15)<<"t_{anal}"<< endl;
    }else{ofile.open(filename, ios::app | ios::out);}

    if (ofile.is_open()){
        ofile << setiosflags(ios::showpoint | ios::uppercase);
        ofile << setw(10) << setprecision(8) << nParticles;
        ofile << setw(15) << setprecision(8) <<t_num/1000.;
        ofile << setw(15) << setprecision(8) << t_anal << "\n";
        ofile.close();
    }else{
        cout << "Error opening file "<<filename << endl;
    }
}

void Sampler::writeStepToFile(int step, int steps){
    // write sampled quantities for each Metropolis step to file

    int i = m_system->getAlphaIndex();
    double alpha = m_system->getWaveFunction()->getParameters()[i];
    int timestep = m_system->getTimeStepIndex();
    double nParticles = m_system->getNumberOfParticles();
    double nDim = m_system->getNumberOfDimensions();
    double nSteps = m_system->getNumberOfMetropolisSteps();
    double dt = m_system->getTimeSteps()[timestep];
    bool calc = m_system->getCalculation();
    bool solv = m_system->getSolver();

    double energy = m_system->getEnergy();

    string type;
    if(calc==true){ type = "numeric"; }
    else if(calc==false){ type = "analytic"; }

    string solver;
    if(solv==true){ solver = "bruteForce"; }
    else if (solv==false){ solver = "importance"; }

    ofstream ofile;
    string filename = "data/c/1c_";
    string arg1 = to_string(int(nParticles));
    string arg2 = to_string(int(nDim));
    string arg3 = to_string(int(nSteps));
    string arg4 = to_string(double(dt));
    string arg5 = to_string(m_system->getWaveFunction()->getParameters()[i]);
    filename.append(solver);
    filename.append("_");
    filename.append(type);
    filename.append("_nParticles_");
    filename.append(arg1);
    filename.append("_nDim_");
    filename.append(arg2);
    filename.append("_nSteps_");
    filename.append(arg3);
    filename.append("_dt_");
    filename.append(arg4);
    filename.append("_alpha_");
    filename.append(arg5);
    filename.append("_.txt");
    if (steps == 0){
        ofile.open(filename, ios::trunc | ios::out);
        ofile << setw(10) << "steps" <<setw(15) << "Energy" << setw(15) << "Variance" << setw(15) << "Error" << endl;
    }
    else{ofile.open(filename, ios::app | ios::out);}

    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(10) << setprecision(8) << step;
    ofile << setw(15) << setprecision(8) << energy;
    ofile << setw(15) << setprecision(8) << m_variance;
    ofile << setw(15) << setprecision(8) << m_error << "\n";
    ofile.close();
}


void Sampler::writeAlphaToFile(){
    // write averaged quantities for each value of alpha to file

    int i = m_system->getAlphaIndex();
    int step = m_system->getMetropolisStep();
    int timestep = m_system->getTimeStepIndex();
    double nParticles = m_system->getNumberOfParticles();
    double nDim = m_system->getNumberOfDimensions();
    double nSteps = m_system->getNumberOfMetropolisSteps();
    double alpha = m_system->getWaveFunction()->getParameters()[i];
    double stepLength = m_system->getStepLength();
    double dt = m_system->getTimeSteps()[timestep];
    bool calc = m_system->getCalculation();
    bool solv = m_system->getSolver();

    cout << "writeAlphaToFile: m_variance = " << m_variance << endl;
    cout << "writeAlphaToFile: m_error = " << m_error << endl;
    cout << "t_num = " << t_num << endl;
    cout << "t_num / CLOCKS_PER_SEC = " << t_num / ((double) CLOCKS_PER_SEC) << endl;

    string type;
    if(calc==true){ type = "numeric"; }
    else if(calc==false){ type = "analytic"; }

    string solver;
    if(solv==true){ solver = "bruteForce"; }
    else if (solv==false){ solver = "importance"; }

    ofstream ofile;
    string filename = "data/c/1c_alpha_";
    string arg1 = to_string(int(nParticles));
    string arg2 = to_string(int(nDim));
    string arg3 = to_string(int(nSteps));
    string arg4 = to_string(stepLength);
    string arg5 = to_string(dt);
    filename.append(solver);
    filename.append("_");
    filename.append(type);
    filename.append("_nPart_");
    filename.append(arg1);
    filename.append("_nDim_");
    filename.append(arg2);
    filename.append("_nSteps_");
    filename.append(arg3);
    filename.append("_stepLength_");
    filename.append(arg4);
    filename.append("_dt_");
    filename.append(arg5);
    filename.append("_.txt");
    if (i == 0){
        ofile.open(filename, ios::trunc | ios::out);
        ofile << setw(10) << "alpha" <<setw(15) << "steps" << setw(15) << "Energy" << setw(15) << "Variance"<< setw(15) << "Error" << setw(15) << "accepted_ratio" << setw(15) << "t_CPU [s]" << "\n"; //<< setw(15) << "Variance" << setw(15) << "Error" << "\n";
    }
    else{ofile.open(filename, ios::app | ios::out);}

    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(10) << setprecision(8) << alpha;
    ofile << setw(15) << setprecision(8) << m_numberOfMetropolisSteps;
    ofile << setw(15) << setprecision(8) << m_cumulativeEnergy / m_stepNumber;
    ofile << setw(15) << setprecision(8) << m_variance;
    ofile << setw(15) << setprecision(8) << m_error;
    ofile << setw(15) << setprecision(8) << m_stepNumber/ ((double) m_numberOfMetropolisSteps);
//    ofile << setw(15) << setprecision(8) << t_num / ((double) m_stepNumber) << "\n";
    ofile << setw(15) << setprecision(8) << t_num  / ((double) CLOCKS_PER_SEC)<< "\n";
    ofile.close();
}

void Sampler::writeTimeStepToFile(){
    // write averaged quantities for each value of time step (dt) to file

    int i = m_system->getTimeStepIndex();
    int step = m_system->getMetropolisStep();
    double nParticles = m_system->getNumberOfParticles();
    double nDim = m_system->getNumberOfDimensions();
    double nSteps = m_system->getNumberOfMetropolisSteps();
    double timestep = m_system->getTimeSteps()[i];
    bool calc = m_system->getCalculation();
    double stepLength = m_system->getStepLength();

    string type;
    if(calc==true){ type = "numeric"; }
    else if(calc==false){ type = "analytic"; }

    ofstream ofile;
    string filename = "data/1c_nParticles_";

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
    filename.append(type);
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

void Sampler::writeOneBodyDensityToFile(){
    // write one-body densities to file

    int i = m_system->getTimeStepIndex();

    int nBins = m_nBins;
    double r_max = m_Rmax;
    bool interacting = m_system->getWaveFunction()->getJastrow(); //adjust for spherical
    int step = m_system->getMetropolisStep();
    double nParticles = m_system->getNumberOfParticles();
    double nDim = m_system->getNumberOfDimensions();
    double nSteps = m_system->getNumberOfMetropolisSteps();
    double timestep = m_system->getTimeSteps()[i];
    bool calc = m_system->getCalculation();
    double stepLength = m_system->getStepLength();

    string type;
    if(calc==true){ type = "numeric"; }
    else if(calc==false){ type = "analytic"; }

    string Jastrow;
    if(interacting==true){ Jastrow = "interacting"; }
    else if(interacting==false){ Jastrow = "non-interacting"; }

    ofstream ofile;
//    string filename = "data/1c_nParticles_";
//    string filename = "data/g/1/final/plusDR_spherical_onebody_dv2_nPart_";
//    string filename = "data/g/1/final/A0.5_elliptical_onebody_dv2_nPart_";
        string filename = "data/g/1/final/SPHERICAL";
    string arg1 = to_string(int(nParticles));
    string arg2 = to_string(int(nDim));
    string arg3 = to_string(int(nSteps));
    string arg4 = to_string(stepLength);
    string arg5 = to_string(m_nBins);
    filename.append(arg1);
    filename.append("_nDim_");
    filename.append(arg2);
    filename.append("_nSteps_");
    filename.append(arg3);
    filename.append("_stepLength_");
    filename.append(arg4);
    filename.append("_nBins_");
    filename.append(arg5);
    filename.append("_");
    filename.append(Jastrow);
    filename.append("_");
    filename.append(type);
    filename.append(".txt");
//    if (i == 0){
        ofile.open(filename, ios::trunc | ios::out);
        ofile << setw(7) << "bin number" << setw(15) << "bin count" <<  "\n";
//    }

    ofile << setiosflags(ios::showpoint | ios::uppercase);

    for(int i = 0; i<nBins; i++){
       ofile << setw(10) << setprecision(8) << (i+1)*(r_max/(double)nBins);
       ofile << setw(15) << setprecision(8) << (double)(m_OneBodyBin[i]/(nSteps)) <<"\n";
    }

    ofile.close();


}
