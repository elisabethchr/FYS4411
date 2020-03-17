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

void Sampler::setNumberOfMetropolisSteps(std::vector<int> numberOfMetropolisSteps) {
    for(int i =0; i<numberOfMetropolisSteps.size(); i++){
          m_numberOfMetropolisSteps.push_back(numberOfMetropolisSteps[i]);
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

    /* Here you should sample all the interesting things you want to measure.
         * Note that there are (way) more than the single one here currently.
         */

    // true = numeric
    bool calc = m_system->getCalculation();

    //calculate cpu-time
    clock_t c_start1 = clock();

    double localEnergy = m_system->getHamiltonian()->
            computeLocalEnergy(m_system->getParticles(), calc);

    clock_t c_end1 = clock();
    t_num += (c_end1-c_start1);

    m_cumulativeEnergy  += localEnergy;
//    m_cumulativeEnergy2 += m_cumulativeEnergy*m_cumulativeEnergy;
    m_cumulativeEnergy2 += localEnergy*localEnergy;


//    clock_t c_start2 = clock();
///*
//    // false = analytic/ not numeric
//    double localEnergyAnalytic = m_system->getHamiltonian()->computeLocalEnergy(m_system->getParticles(), false);
//    m_cumulativeEnergyAnalytic += localEnergyAnalytic;

//    //   cout << t_anal<< "   " << t_num<< endl;
//*/

//    clock_t c_end2 = clock();
//    t_anal += (c_end2-c_start2);

    m_stepNumber++;
}


void Sampler::printOutputToTerminal() {

    int     np = m_system->getNumberOfParticles();
    int     nd = m_system->getNumberOfDimensions();
    int     msIndex = m_system->getnStepsIndex();
    std::vector<int> ms = m_system->getNumberOfMetropolisSteps();
    int     p  = m_system->getWaveFunction()->getNumberOfParameters();
    double  ef = m_system->getEquilibrationFraction();
    std::vector<double> pa = m_system->getWaveFunction()->getParameters();

//    cout << "m_cumulativeenergy = " << m_cumulativeEnergy << endl;
//    cout << "Analytic energy: " << m_cumulativeEnergyAnalytic << endl;

    cout << endl;
    cout << "  -- System info -- " << endl;
    cout << " Number of particles  : " << np << endl;
    cout << " Number of dimensions : " << nd << endl;
//    cout << " Number of Metropolis steps run : 10^" << std::log10(ms) << endl;
//    cout << " Number of equilibration steps  : 10^" << std::log10(std::round(ms*ef)) << endl;
    cout << endl;
    cout << "  -- Wave function parameters -- " << endl;
    cout << " Number of parameters : " << p << endl;
//    for (int i=0; i < p; i++) {
//        cout << " Parameter " << i+1 << " : " << pa.at(i) << endl;
//    }
    cout << endl;
    cout << "  -- Results -- " << endl;
    cout << " Energy : " << m_energy << endl;
    cout << "Variance: " << m_variance << endl;
    cout << "Error: " << m_error << endl;
}

void Sampler::computeAverages() {
    /* Compute the averages of the sampled quantities. You need to think
     * thoroughly through what is written here currently; is this correct?
     */
    //    m_energy = m_cumulativeEnergy / (m_system->getNumberOfMetropolisSteps());
    //    m_energy = m_cumulativeEnergy / (m_system->getNumberOfMetropolisSteps()*m_system->getNumberOfParticles());
    m_energy = m_cumulativeEnergy / ((double) m_stepNumber);
    m_energy2 = m_cumulativeEnergy2 / ((double )m_stepNumber);
    m_energyAnalytic = m_cumulativeEnergyAnalytic / (m_stepNumber);
//    m_variance = m_energy2 - m_energy*m_energy;
    m_variance = (m_cumulativeEnergy2 - m_cumulativeEnergy*m_cumulativeEnergy/(float)m_stepNumber)/((float) m_stepNumber*m_stepNumber);
//    m_variance = (m_cumulativeEnergy2 - m_cumulativeEnergy*m_cumulativeEnergy);
    m_error = pow(abs(m_variance), 0.5);
}


void Sampler::writeTotalToFile(){
    int nParticles = m_system->getNumberOfParticles();
    //    int nDim = m_system->getNumberOfDimensions();
    //    int nSteps = m_system->getNumberOfMetropolisSteps();


    ofstream ofile;
    string filename = "10steps.txt";

    if (nParticles == 1){
        ofile.open(filename, ios::trunc | ios::out);
        ofile << setw(10) << "N_{particles}" <<setw(15) << "t_{num}" << setw(15)<<"t_{anal}"<< endl;
    }else{ofile.open(filename, ios::app | ios::out);}

    //    cout << "m_energy: " << m_energy << endl;
    //    cout << "m_energyAnalytic: " << m_energyAnalytic << endl;
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
    int timestep = m_system->getTimeStepIndex();
    double nParticles = m_system->getNumberOfParticles();
    double nDim = m_system->getNumberOfDimensions();
    int nStepsIndex = m_system->getnStepsIndex();
    std::vector<int> nSteps = m_system->getNumberOfMetropolisSteps();
    double dt = m_system->getTimeSteps()[timestep];
    bool calc = m_system->getCalculation();

    string type;
    if(calc==true){ type = "numeric"; }
    else if(calc==false){ type = "analytic"; }

//    cout << "m_cumulativeEnergy = " << m_cumulativeEnergy << endl;
//    cout << "m_stepNumber = " << m_stepNumber << endl;
//    cout << "m_energy = " << m_cumulativeEnergy / ((double) m_stepNumber) << endl << " " << endl;

    ofstream ofile;
    string filename = "data/c/1c_importance_nParticles_";
//    string filename = "data/test";
    string arg1 = to_string(int(nParticles));
    string arg2 = to_string(int(nDim));
    string arg3 = to_string(int(nSteps[nStepsIndex]));
    string arg4 = to_string(double(dt));
    filename.append(arg1);
    filename.append("_nDim_");
    filename.append(arg2);
    filename.append("_nSteps_");
    filename.append(arg3);
    filename.append("_dt_");
    filename.append(arg4);
    filename.append("_");
    filename.append(type);
    filename.append(".txt");
    if (steps == 0){
        ofile.open(filename, ios::trunc | ios::out);
        ofile << setw(10) << "steps" <<setw(15) << "Energy" << setw(15) << "Variance" << setw(15) << "Error" << endl;
    }
    else{ofile.open(filename, ios::app | ios::out);}

    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(10) << setprecision(8) << step;
    ofile << setw(15) << setprecision(8) << m_cumulativeEnergy / ((double) m_stepNumber);
    ofile << setw(15) << setprecision(8) << m_variance;
    ofile << setw(15) << setprecision(8) << m_error << "\n";
    ofile.close();
}


void Sampler::writeAlphaToFile(){
    int i = m_system->getAlphaIndex();
    int step = m_system->getMetropolisStep();
    double nParticles = m_system->getNumberOfParticles();
    double nDim = m_system->getNumberOfDimensions();
    int nStepsIndex = m_system->getnStepsIndex();
    int nSteps = m_system->getNumberOfMetropolisSteps()[nStepsIndex];
    double alpha = m_system->getWaveFunction()->getParameters()[i];
    bool calc = m_system->getCalculation();
    double stepLength = m_system->getStepLength();

    string type;
    if(calc==true){ type = "numeric"; }
    else if(calc==false){ type = "analytic"; }

    ofstream ofile;
//    string filename = "data/1c_nParticles_";
    string filename = "data/c/brute_alpha_nPart_";
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
        ofile << setw(10) << "alpha" <<setw(15) << "Energy" <<setw(15)<< "Variance" <<
                 setw(15)<<"error"<<"\n";
    }
    else{ofile.open(filename, ios::app | ios::out);}

    //    cout << "m_energy: " << m_energy << endl;
    //    cout << "m_energyAnalytic: " << m_energyAnalytic << endl;

    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(10) << setprecision(8) << alpha;
    ofile << setw(15) << setprecision(8) << m_cumulativeEnergy / m_stepNumber;
    ofile << setw(25) << setprecision(15) << m_variance;
    ofile << setw(25) << setprecision(15) << m_error << "\n";
    ofile.close();
}

void Sampler::writeTimeStepToFile(){
    int i = m_system->getTimeStepIndex();
    int step = m_system->getMetropolisStep();
    double nParticles = m_system->getNumberOfParticles();
    double nDim = m_system->getNumberOfDimensions();
    int nStepsIndex = m_system->getnStepsIndex();
    int nSteps = m_system->getNumberOfMetropolisSteps()[nStepsIndex];
    double timestep = m_system->getTimeSteps()[i];
    bool calc = m_system->getCalculation();
    double stepLength = m_system->getStepLength();

    string type;
    if(calc==true){ type = "numeric"; }
    else if(calc==false){ type = "analytic"; }

    ofstream ofile;
//    string filename = "data/1c_nParticles_";
    string filename = "data/c/timestep_nPart_";
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
        ofile << setw(7) << "dt" << setw(15) << "Energy" << setw(15) << "Variance" << setw(15) << "Error" << setw(15) << "ratio_steps" << "\n";
    }
    else{ofile.open(filename, ios::app | ios::out);}

    //    cout << "m_energy: " << m_energy << endl;
    //    cout << "m_energyAnalytic: " << m_energyAnalytic << endl;

    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(10) << setprecision(8) << timestep;
    ofile << setw(15) << setprecision(8) << m_cumulativeEnergy / ((double) m_stepNumber);
    ofile << setw(15) << setprecision(8) << m_variance;
    ofile << setw(15) << setprecision(8) << m_error;
    ofile << setw(15) << setprecision(8) << m_stepNumber / ((double) m_system->getNumberOfMetropolisSteps()[nStepsIndex]) << "\n";
    ofile.close();

}

void Sampler::writeVarToFile(){
    double acceptanceRatio = m_system->getAcceptanceRatio();
    std::vector<double> alpha = m_system->getWaveFunction()->getParameters();
    int dtIndex = m_system->getTimeStepIndex();
    double dt = m_system->getTimeSteps()[dtIndex];
    double nParticles = m_system->getNumberOfParticles();
    double nDim = m_system->getNumberOfDimensions();
    int nStepsIndex = m_system->getnStepsIndex();
    std::vector<int> nSteps = m_numberOfMetropolisSteps;
//    double numberOfMetropolisSteps = m_system->getMetropolisStep()[nStepsIndex];
    bool calc = m_system->getCalculation();

    string type;
    if(calc==true){ type = "numeric"; }
    else if(calc==false){ type = "analytic"; }

//    cout << "m_cumulativeEnergy = " << m_cumulativeEnergy << endl;
//    cout << "m_stepNumber = " << m_stepNumber << endl;
//    cout << "m_energy = " << m_cumulativeEnergy / ((double) m_stepNumber) << endl << " " << endl;

    ofstream ofile;
    string filename = "data/c/variance_comp/variance_nParticles_";
//    string filename = "data/test";
    string arg1 = to_string(int(nParticles));
    string arg2 = to_string(int(nDim));
    string arg3 = to_string(double(alpha[0]));
    string arg4 = to_string(double(dt));
    filename.append(arg1);
    filename.append("_nDim_");
    filename.append(arg2);
    filename.append("_alpha_");
    filename.append(arg3);
    filename.append("_dt_");
    filename.append(arg4);
    filename.append("_");
    filename.append(type);
    filename.append(".txt");
    if (nStepsIndex == 0){
        ofile.open(filename, ios::trunc | ios::out);
        ofile << setw(10) << "MCsteps" <<setw(15) << "Energy" << setw(15) << "Variance" << setw(15) << "Error"
              <<setw(20) <<"Acceptance rate"<<  endl;
    }
    else{ofile.open(filename, ios::app | ios::out);}

    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(10) << setprecision(8) << nSteps[nStepsIndex];
    ofile << setw(15) << setprecision(8) << m_cumulativeEnergy / ((double) m_stepNumber);
    ofile << setw(15) << setprecision(8) << m_variance;
    ofile << setw(15) << setprecision(8) << m_error;
    ofile << setw(15) << setprecision(8) << acceptanceRatio<< "\n";
    ofile.close();
}
