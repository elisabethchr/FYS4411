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
    double getEnergy()              { return m_energy; }
    double getEnergyDerivative()    { return m_derivativeE; }
    double getEnergyAnalytic()      { return m_energyAnalytic;}
    double getVariance()            { return m_variance;}
    double getError()               { return m_error;}
    void writeTotalToFile();
    void writeStepToFile(int step, int steps);
    void writeAlphaToFile();
    void writeTimeStepToFile();
    void writeOneBodyDensityToFile();



private:
    int     m_numberOfMetropolisSteps = 0;
    int     m_stepNumber = 0;
    double  m_energy = 0;
    double  m_energy2 = 0;
    double  m_energyAnalytic = 0;
    double  m_cumulativeEnergy = 0;
    double  m_cumulativeEnergy2 = 0;
    double  m_cumulativeEnergyAnalytic = 0;
    double  m_variance = 0;
    double  m_error = 0;
    double  m_deltaPsi = 0;
    double  m_deltaEnergy = 0;
    double  m_derivativeE = 0;
    double  m_derivativePsiE = 0;
    double  m_derivativePsi_alpha = 0;
    double  m_cumulativeDeltaPsi = 0;
    double  m_cumulativeDerivativePsiE = 0;


    int                             m_nBins;
    double                          m_Rmax;
    double                          m_Rmin;
    std::vector<double>             m_radii;
    std::vector<double>             m_OneBodyBin;


    clock_t t_num = 0;
    clock_t t_anal = 0;

    class System* m_system = nullptr;
};

//#pragma once
//#include <ctime>

//class Sampler{
//public:
//    Sampler(class System* system);
//    void setNumberOfMetropolisSteps(int steps);
//    void sample(bool acceptedStep);
//    void printOutputToTerminal();
//    void computeAverages();
//    double getEnergy()          { return m_energy; }
//    double getEnergyAnalytic()  { return m_energyAnalytic;}
//    double getVariance()        { return m_variance;}
//    double getError()           { return m_error;}
//    void writeTotalToFile();
//    void writeStepToFile(int step, int steps);
//    void writeAlphaToFile();

//private:
//    int     m_numberOfMetropolisSteps = 0;
//    int     m_stepNumber = 0;
//    double  m_energy = 0;
//    double  m_energy2 = 0;
//    double  m_energyAnalytic = 0;
//    double  m_cumulativeEnergy = 0;
//    double  m_cumulativeEnergy2 = 0;
//    double  m_cumulativeEnergyAnalytic = 0;
//    double  m_variance = 0;
//    double  m_error = 0;

//    clock_t t_num = 0;
//    clock_t t_anal = 0;

//    class System* m_system = nullptr;
//};
