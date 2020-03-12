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

    clock_t c_start1 = clock();

    double localEnergy = m_system->getHamiltonian()->
            computeLocalEnergy(m_system->getParticles());
    clock_t c_end1 = clock();
    t_num += (c_end1-c_start1);
    m_cumulativeEnergy  += localEnergy;

    //cout << "c_start1  " << c_start1 << "c_end1  " << c_end1 << endl;



    clock_t c_start2 = clock();

    double localEnergyAnalytic = m_system->getHamiltonian()->computeLocalEnergy_analytic(m_system->getParticles());
    clock_t c_end2 = clock();
    t_anal += (c_end2-c_start2);
    m_cumulativeEnergyAnalytic += localEnergyAnalytic;

//   cout << t_anal<< "   " << t_num<< endl;



    m_stepNumber++;
}

void Sampler::printOutputToTerminal() {
    int     np = m_system->getNumberOfParticles();
    int     nd = m_system->getNumberOfDimensions();
    int     ms = m_system->getNumberOfMetropolisSteps();
    int     p  = m_system->getWaveFunction()->getNumberOfParameters();
    double  ef = m_system->getEquilibrationFraction();
    std::vector<double> pa = m_system->getWaveFunction()->getParameters();

    cout << "m_cumulativeenergy = " << m_cumulativeEnergy << endl;
    cout << "Analytic energy: " << m_cumulativeEnergyAnalytic << endl;

    cout << endl;
    cout << "  -- System info -- " << endl;
    cout << " Number of particles  : " << np << endl;
    cout << " Number of dimensions : " << nd << endl;
    cout << " Number of Metropolis steps run : 10^" << std::log10(ms) << endl;
    cout << " Number of equilibration steps  : 10^" << std::log10(std::round(ms*ef)) << endl;
    cout << endl;
    cout << "  -- Wave function parameters -- " << endl;
    cout << " Number of parameters : " << p << endl;
    for (int i=0; i < p; i++) {
        cout << " Parameter " << i+1 << " : " << pa.at(i) << endl;
    }
    cout << endl;
    cout << "  -- Reults -- " << endl;
    cout << " Energy : " << m_energy << endl;
    cout << " Analytic energy: " << m_energyAnalytic << endl;
    cout << "Numerical CPU time[ms]  "<< t_num << endl;
    cout << "Analytical CPU time[ms]  " << t_anal << endl;
}

void Sampler::computeAverages() {
    /* Compute the averages of the sampled quantities. You need to think
     * thoroughly through what is written here currently; is this correct?
     */
    //    m_energy = m_cumulativeEnergy / (m_system->getNumberOfMetropolisSteps());
    //    m_energy = m_cumulativeEnergy / (m_system->getNumberOfMetropolisSteps()*m_system->getNumberOfParticles());

    m_energy = m_cumulativeEnergy / m_stepNumber;
    m_energyAnalytic = m_cumulativeEnergyAnalytic / (m_stepNumber);
}

void Sampler::writeToFile(int step, int steps){
    double nParticles = m_system->getNumberOfParticles();
    double nDim = m_system->getNumberOfDimensions();
    double nSteps = m_system->getNumberOfMetropolisSteps();


    ofstream ofile;
    string filename = "b_nParticles_";
    string arg1 = to_string(int(nParticles));
    string arg2 = to_string(int(nDim));
    string arg3 = to_string(int(nSteps));
    filename.append(arg1);
    filename.append("_nDim");
    filename.append(arg2);
    filename.append("_nSteps_");
    filename.append(arg3);
    filename.append(".txt");
    if (steps == 0){ofile.open(filename, ios::trunc | ios::out);}
    else{ofile.open(filename, ios::app | ios::out);}

//    cout << "m_energy: " << m_energy << endl;
//    cout << "m_energyAnalytic: " << m_energyAnalytic << endl;
    if (ofile.is_open()){
        ofile << setiosflags(ios::showpoint | ios::uppercase);
        ofile << setw(10) << setprecision(8) << step;
        ofile << setw(15) << setprecision(8) << m_cumulativeEnergy / m_stepNumber;
        ofile << setw(15) << setprecision(8) << m_cumulativeEnergyAnalytic / m_stepNumber << "\n";
        ofile.close();
    }else{
        cout << "Error opening file "<<filename << endl;
    }

}