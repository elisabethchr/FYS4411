#pragma once
#include <vector>
#include <armadillo>

class System {
public:
//    bool metropolisStep             (unsigned gen);
    bool bruteForce                 ();     // Standard Metropolis
    bool importanceSampling         ();     // Metropolis-Hastings
    void runMetropolisSteps         (int RBM_cycles, std::vector<int> numberOfMetropolisSteps);
    void setNumberHiddenNodes       (int n_hidden){ m_numberHiddenNodes = n_hidden; }
    void setNumberVisibleNodes      (int n_visible){ m_numberVisibleNodes = n_visible; }
    void setNumberParticles         (int nPart);
    void setNumberDimensions        (int nDim);
    void setStepLength              (double stepLength);
    void setEquilibrationFraction   (double equilibrationFraction);
    void setWaveFunctionValue       (double waveFunction);
    void setTimeSteps               (arma::vec timesteps){ m_timesteps = timesteps; }
    void setSolver                  (bool bruteForce){ m_solver = bruteForce; }
    void setLearningRate            (double eta){ m_eta = eta; }
    void setHamiltonian             (class Hamiltonian* hamiltonian);
    void setWaveFunction            (class WaveFunction* waveFunction);
    void setOptimizer               (class Optimizer* optimizer);
    void setInitialState            (class InitialState* initialState);
//    void setSampler                 (class Sampler* sampler);
    class WaveFunction*             getWaveFunction()   { return m_waveFunction; }
    class Hamiltonian*              getHamiltonian()    { return m_hamiltonian; }
    class Sampler*                  getSampler()        { return m_sampler; }
    class Optimizer*                getOptimizer()      { return m_optimizer; }
    std::vector<class Particle*>    getParticles()      { return m_particles; }
    int getNumberHiddenNodes()          { return m_numberHiddenNodes; }
    int getNumberVisibleNodes()         { return m_numberVisibleNodes; }
    int getNumberParticles()            { return m_numberParticles; }
    int getNumberDimensions()           { return m_numberDimensions; }
    int getNumberOfMetropolisSteps()    { return m_numberOfMetropolisSteps; }
    int getNumberRBMcycles()            { return m_RBMcycles; }
    int getTimeStepIndex()              { return m_timestep; }
    int getRBMstep()                    { return m_RBMstep; }
    int getMetropolisStep()             { return m_MCstep; }
    int getSampleStep()                 { return m_sampleStep; }
    bool getSolver()                    { return m_solver; }
    double getEquilibrationFraction()   { return m_equilibrationFraction; }
    double getStepLength()              { return m_stepLength; }
    double getLearningRate()            { return m_eta; }
    double getEnergyDerivative()        { return m_derivativeE; }
    double getWaveFunctionValue()       { return m_wfValue; }
    double getUniform(double min, double max)    { std::uniform_real_distribution<float> gen(min, max); return gen(m_randomengine); }
    double getGaussian(double mean, double std)  { std::normal_distribution<float> gen(mean, std); return gen(m_randomengine); }
    arma::vec getTimeSteps()           { return m_timesteps; }


private:
    int                             m_numberHiddenNodes = 0;
    int                             m_numberVisibleNodes = 0;
    int                             m_numberOfMetropolisSteps = 0;
    int                             m_RBMcycles = 0;
    int                             m_numberParticles = 0;
    int                             m_numberDimensions = 0;
    int                             m_timestep = 0;
    int                             m_acceptedSteps = 0;
    int                             m_RBMstep = 0;
    int                             m_MCstep;
    int                             m_sampleStep;
    double                          m_acceptedSteps_ratio = 0.0;
    double                          m_equilibrationFraction = 0.0;
    double                          m_stepLength = 0.1;
    double                          m_derivativeE = 0.0;
    double                          m_wfValue = 0.0;
    double                          m_eta = 0.0;        // learning rate
    bool                            m_solver;
    class WaveFunction*             m_waveFunction = nullptr;
    class Hamiltonian*              m_hamiltonian = nullptr;
    class InitialState*             m_initialState = nullptr;
    class Sampler*                  m_sampler = nullptr;
    class Optimizer*                m_optimizer = nullptr;
    std::vector<class Particle*>    m_particles = std::vector<class Particle*>();
    arma::vec                       m_timesteps;

    std::mt19937_64 m_randomengine;

};
