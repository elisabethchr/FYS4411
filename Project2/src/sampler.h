#pragma once
#include <ctime>
#include <armadillo>
#include <iostream>

using namespace std;

class Sampler{
public:
    Sampler(class System* system);
    void setNumberOfMetropolisSteps(int steps);
    void sample(bool acceptedStep);
    void printOutputToTerminal();
    void computeAverages();
    void optimizeWeights();
//    double getEnergy()              { return m_energy; }
    void writeToFile(int step, int steps, string filename);
    void writeStepToFile(int step, int steps, string filename);
    void writeTimeStepToFile(string filename);



private:
    int     m_numberOfMetropolisSteps = 0;
    int     m_nh, m_nv;
    double  m_stepNumber = 0;
    double  m_energy;   // mean energy
    double  m_energy2;  // mean energy^2
    double  m_energyAnalytic = 0;
    double  m_cumulativeEnergy = 0;
    double  m_cumulativeEnergy2 = 0;
    double  m_cumulativeEnergyAnalytic = 0;
    double  m_variance = 0;
    double  m_error = 0;
    double  m_deltaPsi = 0;
    double  m_deltaEnergy = 0;  // change in energy
    double  m_deltaVariance = 0; // variance in energy for a certain Metropolis step
    double  m_derivativeE = 0;  // derivative of energy wrt. variational parameter alpha
    arma::vec  m_gradE;
    arma::vec  m_EdPsi;   // (derivative of psi) * (change in energy)
    arma::vec  m_dPsi;  // derivative of psi wrt. variational parameter alpha
    arma::vec  m_cumulative_dPsi;
    arma::vec  m_cumulative_EdPsi;

    clock_t t_num = 0;
    clock_t t_anal = 0;

    class System* m_system = nullptr;
};
