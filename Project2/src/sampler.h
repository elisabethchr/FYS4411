#pragma once
#include <ctime>
#include <armadillo>

class Sampler{
public:
    Sampler(class System* system);
    void setNumberOfMetropolisSteps(int steps);
    void setOneBodyDensity          (int nBins, double r_max, double r_min);
    void sample(bool acceptedStep);
    void sampleOneBodyDensity();
    void printOutputToTerminal();
    void computeAverages();
    void optimizeWeights();
    double getEnergy()              { return m_energy; }
    double getEnergyDerivative()    { return m_derivativeE; }
    double getEnergyAnalytic()      { return m_energyAnalytic;}
    double getVariance()            { return m_variance;}
    double getError()               { return m_error;}
    void writeToFile();
    void writeStepToFile(int step, int steps);
    void writeTimeStepToFile();



private:
    int     m_numberOfMetropolisSteps = 0;
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
    double  m_derivativeE = 0;  // derivative of energy wrt. variational parameter alpha
    arma::vec  m_gradE;
    arma::vec  m_EdPsi = 0;   // (derivative of psi) * (change in energy)
    arma::vec  m_dPsi = 0;  // derivative of psi wrt. variational parameter alpha
    arma::vec  m_cumulative_dPsi = 0;
    arma::vec  m_cumulative_EdPsi = 0;

    int                             m_nBins;
    double                          m_Rmax;
    double                          m_Rmin;
    std::vector<double>             m_radii;
    std::vector<double>             m_OneBodyBin;


    clock_t t_num = 0;
    clock_t t_anal = 0;

    class System* m_system = nullptr;
};
