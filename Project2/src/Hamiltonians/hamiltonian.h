#pragma once
#include <vector>
#include <armadillo>

using namespace arma;

class Hamiltonian {
public:
    Hamiltonian(class System* system);

    virtual double getOmega() = 0;
    virtual double computeLocalEnergy(arma::vec X) = 0;

    virtual vec computeQuantumForce(arma::vec X, int i) = 0;
    virtual vec computeGradientPsi(arma::vec X) = 0;


protected:
    class System* m_system = nullptr;
};
