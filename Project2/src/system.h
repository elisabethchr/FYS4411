#pragma once
#include <vector>
#include <armadillo>

class System {
public:
//    bool metropolisStep             (unsigned gen);
    bool metropolisStep             ();
    bool importanceSampling         ();     // Metropolis-Hastings
    void runMetropolisSteps         (std::vector<int> numberOfMetropolisSteps);
    void setNumberHiddenNodes       (int n_hidden);
    void setNumberVisibleNodes      (int n_visible);
    void setStepLength              (double stepLength);
    void setEquilibrationFraction   (double equilibrationFraction);
    void setWaveFunctionValue       (double waveFunction);
    void setTimeSteps               (std::vector<double> timesteps){ m_timesteps = timesteps; }
    void setSolver                  (bool bruteForce){ m_solver = bruteForce; }
    void setHamiltonian             (class Hamiltonian* hamiltonian);
    void setWaveFunction            (class WaveFunction* waveFunction);
    void setInitialState            (class InitialState* initialState);
    void setEnergy                  (double localEnergy){m_energy = localEnergy;}
    class WaveFunction*             getWaveFunction()   { return m_waveFunction; }
    class Hamiltonian*              getHamiltonian()    { return m_hamiltonian; }
    class Sampler*                  getSampler()        { return m_sampler; }
    std::vector<class Particle*>    getParticles()      { return m_particles; }
    int getNumberHiddenNodes()          { return m_numberHiddenNodes; }
    int getNumberVisibleNodes()         { return m_numberVisibleNodes; }
    int getNumberOfMetropolisSteps()    { return m_numberOfMetropolisSteps; }
    int getAlphaIndex()                 { return m_alpha; }
    int getTimeStepIndex()              { return m_timestep; }
    int getMetropolisStep()             { return m_stepMetropolis; }
    bool getSolver()                    { return m_solver; }
    double getEquilibrationFraction()   { return m_equilibrationFraction; }
    double getStepLength()              { return m_stepLength; }
    double getEnergy()                  { return m_energy; }
    double getEnergyDerivative()        { return m_derivativeE; }
    double getWaveFunctionValue()       { return m_wfValue; }
    double getUniform(double min, double max)    { std::uniform_real_distribution<float> gen(min, max); return gen(m_seed); }
    double getGaussian(double mean, double std)  { std::normal_distribution<float> gen(mean, std); return gen(m_seed); }
    std::vector<double> getTimeSteps()           { return m_timesteps; }


private:
    int                             m_numberHiddenNodes = 0;
    int                             m_numberVisibleNodes = 0;
    int                             m_numberOfMetropolisSteps = 0;
    int                             m_alpha = 0;
    int                             m_timestep = 0;
    int                             m_acceptedSteps = 0;
    int                             m_stepMetropolis = 0;
    int                             m_stepImportance = 0;
    int                             m_MCstep = 0;
    double                          m_acceptedSteps_ratio = 0.0;
    double                          m_equilibrationFraction = 0.0;
    double                          m_stepLength = 0.1;
    double                          m_energy = 0.0;
    double                          m_derivativeE = 0.0;
    double                          m_wfValue = 0.0;
    bool                            m_solver;
    class WaveFunction*             m_waveFunction = nullptr;
    class Hamiltonian*              m_hamiltonian = nullptr;
    class InitialState*             m_initialState = nullptr;
    class Sampler*                  m_sampler = nullptr;
    std::vector<class Particle*>    m_particles = std::vector<class Particle*>();
    std::vector<double>             m_timesteps;

    std::mt19937_64 m_seed;

};
