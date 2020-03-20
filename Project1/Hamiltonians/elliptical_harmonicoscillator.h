#pragma once
#include "hamiltonian.h"
#include <vector>

class EllipticalHarmonicOscillator : public Hamiltonian {
public:
    EllipticalHarmonicOscillator(System* system, double omega);
    double computeLocalEnergy(std::vector<Particle*> particles, bool type);
    vec computeQuantumForce(std::vector<Particle*> particles, int i);

private:
    double m_omega = 0;
};
