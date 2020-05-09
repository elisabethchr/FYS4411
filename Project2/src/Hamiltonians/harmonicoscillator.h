#pragma once
#include "hamiltonian.h"
#include <vector>

class HarmonicOscillator : public Hamiltonian {
public:
    HarmonicOscillator(System* system, double omega, bool interaction);

    double getOmega() { return m_omega; }
    double computeLocalEnergy(arma::vec X);     // compute local energy as a function of visible nodes X

    vec computeQuantumForce(arma::vec X, int i);
    vec computeGradientPsi(arma::vec X){}

private:
    double m_omega = 0;     // oscillator frequency
    bool m_interaction;     // boolean for whether or not to include interaction
};

