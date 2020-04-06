#pragma once
#include <vector>
#include <armadillo>

class WaveFunction {
public:
    WaveFunction(class System* system);
    int     getNumberOfParameters()     { return m_numberOfParameters; }
    int     getNumberOfParametersBeta() { return m_numberOfParametersBeta; }
    double  getHardCoreDiameter()       { return m_hard_core_diameter; }
    double  getGamma()                  { return m_gamma; }
    bool getJastrow()                   {return m_Jastrow;}
    std::vector<double> getParameters()     { return m_parameters; }
    std::vector<double> getParametersBeta() { return m_beta; }

    virtual double evaluate(std::vector<class Particle*> particles) = 0;        //"virtual" implements function elsewhere, "virtual = 0" forces to implemented elsewhere
    virtual arma::vec computeGradient(std::vector<class Particle*> particles, int i) = 0;
    virtual double computeLaplacian(std::vector<class Particle*> particles) = 0;
    virtual double computeDoubleDerivative_numeric(std::vector<class Particle*> particles) = 0;
    virtual double computeDoubleDerivative_analytic(std::vector<class Particle*> particles) = 0;
    virtual double getDistance(std::vector<class Particle *> particles, int particle_i, int particle_j) = 0;
    virtual double computeDerivativePsi_alpha(std::vector<class Particle*> particles) = 0;

protected:
    int     m_numberOfParameters = 0;
    int     m_numberOfParametersBeta = 0;
    double  m_hard_core_diameter = 0;
//    double  m_gamma = 2.83843;
    double  m_gamma = 1.0;
    bool m_Jastrow;

    std::vector<double> m_parameters = std::vector<double>();
    std::vector<double> m_beta = std::vector<double>();
    class System* m_system = nullptr;

    virtual double SingleParticleFunction(std::vector<class Particle *> particles, int particle) = 0;
    virtual double correlationWaveFunction(std::vector<class Particle* > particles, double distance) = 0;
    virtual double computeLaplacian_SPF(int particle_i) = 0;
    virtual arma::vec computeGradient_SPF(int particle_i) = 0;
    virtual double computeLaplacian_correlation(std::vector<class Particle *> particles, int particle_i) = 0;
    virtual arma::vec computeGradient_correlation(std::vector<class Particle *> particles, int particle_i) = 0;
    virtual double computeDerivative_correlation(int particle_i, int particle_j) = 0;
    virtual double computeDoubleDerivative_correlation(int particle_i, int particle_j) = 0;

};
