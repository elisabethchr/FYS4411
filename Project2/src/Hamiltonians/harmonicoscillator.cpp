#include "harmonicoscillator.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"
#include <vector>
#include <armadillo>

using std::cout;
using std::endl;
using namespace arma;

HarmonicOscillator::HarmonicOscillator(System* system, double omega) :
    Hamiltonian(system) {
    assert(omega > 0);
    m_omega  = omega;
}

double HarmonicOscillator::computeLocalEnergy(std::vector<Particle*> particles, bool type) {
    // function for computing the local energy of a spherical system (no perturbation)

    return 0;
}

vec HarmonicOscillator::computeQuantumForce(std::vector<Particle *> particles, int i){
    // Compute the drift force experienced by particles


//    int nDim = m_system->getNumberOfDimensions();
    int nDim = 2;
    vec force(nDim); force.zeros();
    vec gradient = m_system->getWaveFunction()->computeGradient(particles, i);

    // obtain current value for the wavefunction
    double psi = m_system->getWaveFunctionValue();

    force = 2*gradient / psi;

    return force;
}
