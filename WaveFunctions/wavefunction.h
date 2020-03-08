#pragma once
#include <vector>


class WaveFunction {
public:
    WaveFunction(class System* system);
    int     getNumberOfParameters() { return m_numberOfParameters; }
    std::vector<double> getParameters() { return m_parameters; }
    virtual double evaluate(std::vector<class Particle*> particles) = 0;        //"virtual" implements function elsewhere, "virtual = 0" forces to implemented elsewhere
    virtual double computeDoubleDerivative(std::vector<class Particle*> particles) = 0;

protected: //remove later: Accessible to derived classes, but not public
    int     m_numberOfParameters = 0;
    std::vector<double> m_parameters = std::vector<double>(); //remove later: Parameters for the specific wavefunction in question. e.g. for the derived class simplegaussian, it's where the alpha's are stored
    class System* m_system = nullptr;
};

