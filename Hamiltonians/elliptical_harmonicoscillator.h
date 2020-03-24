#pragma once
#include "hamiltonian.h"
#include <vector>

class EllipticalHarmonicOscillator : public Hamiltonian {
public:
    EllipticalHarmonicOscillator(System* system, double omega);
    double computeLocalEnergy(std::vector<Particle*> particles, bool type);
    vec computeQuantumForce(std::vector<Particle*> particles, int i);       //must implement function which is written in project description!
    // -> (computeQuantumForce function might be better in WaveFunctions along )

private:
    double m_omega = 0;
};
