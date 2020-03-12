#pragma once
#include <vector>
#include <armadillo>

class System {
public:
//    bool metropolisStep             (arma::mat dx_mat);
    bool metropolisStep             ();
    bool importanceSampling         ();     // Metropolis-Hastings
    void runMetropolisSteps         (int numberOfMetropolisSteps);
    void setNumberOfParticles       (int numberOfParticles);
    void setNumberOfDimensions      (int numberOfDimensions);
    void setStepLength              (double stepLength);
    void setEquilibrationFraction   (double equilibrationFraction);
    void setCalculation             (bool numeric){m_numeric=numeric; }
    void setHamiltonian             (class Hamiltonian* hamiltonian);
    void setWaveFunction            (class WaveFunction* waveFunction);
    void setInitialState            (class InitialState* initialState);
    class WaveFunction*             getWaveFunction()   { return m_waveFunction; }
    class Hamiltonian*              getHamiltonian()    { return m_hamiltonian; }
    class Sampler*                  getSampler()        { return m_sampler; }
    std::vector<class Particle*>    getParticles()      { return m_particles; }
    int getNumberOfParticles()          { return m_numberOfParticles; }
    int getNumberOfDimensions()         { return m_numberOfDimensions; }
    int getNumberOfMetropolisSteps()    { return m_numberOfMetropolisSteps; }
    int getAlphaIndex()                 { return m_alpha; }
    int getMetropolisStep()             { return m_stepMetropolis; }
    double getEquilibrationFraction()   { return m_equilibrationFraction; }
    double getStepLength()              { return m_stepLength; }
    double getEnergy()                  { return m_energy; }
    bool getCalculation()               { return m_numeric; }


private:
    int                             m_numberOfParticles = 0;
    int                             m_numberOfDimensions = 0;
    int                             m_numberOfMetropolisSteps = 0;
    int                             m_alpha = 0;
    int                             m_stepMetropolis = 0;
    int                             m_stepImportance = 0;
    double                          m_equilibrationFraction = 0.0;
//    double                          m_stepLength = 0.1;
    double                          m_stepLength = 0.1;
    double                          m_energy = 0.0;
    bool                            m_numeric;
    class WaveFunction*             m_waveFunction = nullptr;
    class Hamiltonian*              m_hamiltonian = nullptr;
    class InitialState*             m_initialState = nullptr;
    class Sampler*                  m_sampler = nullptr;
    std::vector<class Particle*>    m_particles = std::vector<class Particle*>();
};

