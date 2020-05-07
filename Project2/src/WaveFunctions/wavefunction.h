#pragma once
#include <vector>
#include <armadillo>

class WaveFunction {
public:
    WaveFunction(class System* system);
    virtual double evaluate(std::vector<class Particle*> particles) = 0;        //"virtual" implements function elsewhere, "virtual = 0" forces to implemented elsewhere
    virtual arma::vec computeGradient(std::vector<class Particle*> particles, int i) = 0;
    virtual double computeLaplacian(std::vector<class Particle*> particles) = 0;
    virtual double computeDoubleDerivative_numeric(std::vector<class Particle*> particles) = 0;
    virtual double getDistance(std::vector<class Particle *> particles, int particle_i, int particle_j) = 0;
    virtual double computeDerivativePsi_alpha(std::vector<class Particle*> particles) = 0;

protected:
    class System* m_system = nullptr;

};
