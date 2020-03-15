#pragma once
#include <vector>
#include <armadillo>

class WaveFunction {
public:
    WaveFunction(class System* system);
    int     getNumberOfParameters() { return m_numberOfParameters; }
    std::vector<double> getParameters() { return m_parameters; }
    virtual double evaluate(std::vector<class Particle*> particles) = 0;        //"virtual" implements function elsewhere, "virtual = 0" forces to implemented elsewhere
    virtual arma::vec computeGradient(std::vector<class Particle*> particles, int i) = 0;
    virtual double computeDoubleDerivative_numeric(std::vector<class Particle*> particles) = 0;
    virtual double computeDoubleDerivative_analytic(std::vector<class Particle*> particles) = 0;

protected:
    int     m_numberOfParameters = 0;
    std::vector<double> m_parameters = std::vector<double>();
    class System* m_system = nullptr;
};

