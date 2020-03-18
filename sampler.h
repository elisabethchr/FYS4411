#pragma once
#include <ctime>
#include <vector>

class Sampler{
public:
    Sampler(class System* system);
    void setNumberOfMetropolisSteps(std::vector<int> numberOfMetropolisSteps);
    void sample(bool acceptedStep);
    void printOutputToTerminal();
    void computeAverages();
    double getEnergy()          { return m_energy; }
    double getEnergyAnalytic()  { return m_energyAnalytic;}
    double getVariance()        { return m_variance;}
    double getError()           { return m_error;}
    void writeTotalToFile();
    void writeStepToFile(int step, int steps);
    void writeAlphaToFile();
    void writeTimeStepToFile();
    void writeVarToFile(bool bruteForce);
    void write3DToFile(bool bruteForce);

private:
    std::vector<int>     m_numberOfMetropolisSteps;
    int     m_stepNumber = 0;
    double  m_energy = 0;
    double  m_energy2 = 0;
    double  m_energyAnalytic = 0;
    double  m_cumulativeEnergy = 0;
    double  m_cumulativeEnergy2 = 0;
    double  m_cumulativeEnergyAnalytic = 0;
    double  m_variance = 0;
    double  m_error = 0;

    clock_t t_num = 0;
    clock_t t_anal = 0;

    class System* m_system = nullptr;
};
