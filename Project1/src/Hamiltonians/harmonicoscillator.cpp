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

    double potentialEnergy = 0;
    double kineticEnergy   = 0;
    int dim = m_system->getNumberOfDimensions();
    int nPart= m_system->getNumberOfParticles();

    // compute local energy
    for (int i = 0; i<nPart;i++){
        for(int j = 0; j<dim; j++){
            double coordinate = particles.at(i)->getPosition().at(j);
            potentialEnergy+= 0.5*m_omega*m_omega*coordinate*coordinate;
        }
    }

    // obtain current value for the wavefunction
    double psi = m_system->getWaveFunctionValue();

    // compute kinetic energy
    if(type == true){ kineticEnergy -= (1/psi)*0.5*m_system->getWaveFunction()->computeDoubleDerivative_numeric(particles); }
    else if(type == false){ kineticEnergy = (1/psi)*0.5*m_system->getWaveFunction()->computeDoubleDerivative_analytic(particles); }

    return potentialEnergy + kineticEnergy;
}

vec HarmonicOscillator::computeQuantumForce(std::vector<Particle *> particles, int i){
    // Compute the drift force experienced by particles

    int nDim = m_system->getNumberOfDimensions();
    vec force(nDim); force.zeros();
    vec gradient = m_system->getWaveFunction()->computeGradient(particles, i);

    // obtain current value for the wavefunction
    double psi = m_system->getWaveFunctionValue();


    force = 2*gradient / psi;


    return force;
}
