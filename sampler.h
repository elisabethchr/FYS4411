#pragma once

class Sampler {
public:
    Sampler(class System* system);
    void setNumberOfMetropolisSteps(int steps);
    void sample(bool acceptedStep);
    void printOutputToTerminal();
    void computeAverages();
    double getEnergy()          { return m_energy; }
    void writeStepToFile(int step, int steps);
    void writeTotalToFile();

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
    clock_t t_num = 0;
    clock_t t_anal = 0;

    class System* m_system = nullptr;
};
