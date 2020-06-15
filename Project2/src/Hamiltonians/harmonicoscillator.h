#pragma once
#include "../system.h"
#include "hamiltonian.h"
#include <vector>

class HarmonicOscillator : public Hamiltonian {
public:
    HarmonicOscillator(System* system, double omega, bool interaction);

    double getOmega() { return m_omega; }
    double computeLocalEnergy(arma::vec X);     // compute local energy as a function of visible nodes X
    double interaction(arma::vec X, int nx, int dim){};

    vec computeLocalEnergyGradient();
    vec computeQuantumForce(arma::vec X, int i);
    vec computeGradientPsi(arma::vec X){}

private:
    double m_omega = 0;     // oscillator frequency
    double m_sigma;         // error of gaussian distribution
    double m_sigma2;        // variance of gaussian distribution
    bool m_interaction;     // boolean for whether or not to include interaction
    int m_nv;
    int m_nh;
};

