#pragma once
#include "hamiltonian.h"
#include <vector>

// class <name1> : public <name2> --> class name1 inherits all methods and member variables from class name2.
// In this case since the class Hamiltonian is defined with System as argument (with getWavefunction, getHamiltonian, getSampler and getParticles as function calls),
// then class HarmonicOscillator can now also get access to the functions defined within System through Hamiltonian.
// And from class System, we have access to all other member variables and functions of specific classes.
// E.g. to get position of particles when in class HarmonicOscillator: std::vector<particle *> particles = m_particles[i].getPosition()[0] (m_particles is explicitly defined in System).
// To get access to parameters (which is in the WaveFunction class): double alpha = m_system->getWaveFunction()->getParameters[0]. (In getWaveFunction we now gain access to all
// member functions/calls and member variables within the WaveFunction class).
class HarmonicOscillator : public Hamiltonian {
public:
    HarmonicOscillator(System* system, double omega);
    double computeLocalEnergy(std::vector<Particle*> particles);
    double computeLocalEnergy_analytic(std::vector<Particle*> particles);

private:
    double m_omega = 0;
};

