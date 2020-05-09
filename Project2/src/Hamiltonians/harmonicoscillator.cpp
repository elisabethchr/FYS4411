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

HarmonicOscillator::HarmonicOscillator(System* system, double omega, bool interaction) :
    Hamiltonian(system) {
    assert(omega > 0);
    m_omega  = omega;
    m_interaction = interaction;
}

/* Function for computing the local energy of a spherical system (no perturbation) */
double HarmonicOscillator::computeLocalEnergy(arma::vec X) {


    return 0;
}

/* Compute the drift force experienced by particles */
vec HarmonicOscillator::computeQuantumForce(arma::vec X, int i){

/*
//    int nDim = m_system->getNumberOfDimensions();
    int nDim = 2;
    vec force(nDim); force.zeros();
    vec gradient = m_system->getWaveFunction()->computeGradient(particles, i);

    // obtain current value for the wavefunction
    double psi = m_system->getWaveFunctionValue();

    force = 2*gradient / psi;

    return force;
*/
    arma::vec a = zeros(2);
    return a;
}
