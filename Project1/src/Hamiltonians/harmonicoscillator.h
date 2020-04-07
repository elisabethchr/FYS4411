#pragma once
#include "hamiltonian.h"
#include <vector>

class HarmonicOscillator : public Hamiltonian {
public:
    HarmonicOscillator(System* system, double omega);
    double computeLocalEnergy(std::vector<Particle*> particles, bool type);
    double getOmega() { return m_omega; }
    vec computeQuantumForce(std::vector<Particle*> particles, int i);
    vec computeGradientPsi(std::vector<Particle*> particles){}

private:
    double m_omega = 0;     // oscillator frequency
};

