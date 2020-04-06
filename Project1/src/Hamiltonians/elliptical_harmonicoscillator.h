#pragma once
#include "hamiltonian.h"
#include <vector>

class EllipticalHarmonicOscillator : public Hamiltonian {
public:
    EllipticalHarmonicOscillator(System* system, double omega);
    double computeLocalEnergy(std::vector<Particle*> particles, bool type);
    double getOmega(){ return m_omega; }
    vec computeQuantumForce(std::vector<Particle*> particles, int i);       //must implement function which is written in project description!
    vec computeGradientPsi(std::vector<Particle*> particles);

private:
    double m_omega = 0;
};
