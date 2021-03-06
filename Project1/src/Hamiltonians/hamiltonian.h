#pragma once
#include <vector>
#include <armadillo>

using namespace arma;

class Hamiltonian {
public:
    Hamiltonian(class System* system);
    virtual double computeLocalEnergy(std::vector<class Particle*> particles, bool type) = 0;
    virtual double getOmega() = 0;
    virtual vec computeQuantumForce(std::vector<class Particle*> particles, int i) = 0;
    virtual vec computeGradientPsi(std::vector<Particle*> particles) = 0;


protected:
    class System* m_system = nullptr;
};
