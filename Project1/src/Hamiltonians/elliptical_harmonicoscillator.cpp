#include "elliptical_harmonicoscillator.h"
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

EllipticalHarmonicOscillator::EllipticalHarmonicOscillator(System* system, double omega) :
    Hamiltonian(system) {
    assert(omega > 0);
    m_omega  = omega;
}

double EllipticalHarmonicOscillator::computeLocalEnergy(std::vector<Particle*> particles, bool type) {
    // function for computing the local energy of an elliptical system (with perturbation)

    double coordinate;
    double potentialEnergy = 0.0;
    double kineticEnergy   = 0.0;
    double beta = m_system->getWaveFunction()->getParametersBeta()[0];
    double gamma = m_system->getWaveFunction()->getGamma();
    double hbar = 1.0;

    int nDim = m_system->getNumberOfDimensions();
    int nPart= m_system->getNumberOfParticles();

    std::vector<double> pos;

    // compute potential energy
    for (int i = 0; i<nPart;i++){
        pos = particles[i]->getPosition();
        for(int j = 0; j<nDim; j++){

            // add perturbation to z-direction
            if (j==2){ coordinate = gamma*pos[j]; }

            else { coordinate = pos[j]; }

            potentialEnergy+= 0.5*hbar*m_omega*m_omega*coordinate*coordinate;
        }
    }

    // obtain current value for the wavefunction
    double psi = m_system->getWaveFunctionValue();

    // compute kinetic energy
    if(type == true){ kineticEnergy -= (1/psi)*0.5*(hbar*hbar/M)*m_system->getWaveFunction()->computeDoubleDerivative_numeric(particles); }
    else if(type == false){ kineticEnergy -= 0.5*(hbar*hbar/M)*m_system->getWaveFunction()->computeDoubleDerivative_analytic(particles); }


    return potentialEnergy + kineticEnergy;
}

vec EllipticalHarmonicOscillator::computeQuantumForce(std::vector<Particle *> particles, int i){
     // compute the drift force experienced by particles

    int nDim = m_system->getNumberOfDimensions();
    double psi = m_system->getWaveFunctionValue();

    vec force(nDim); force.zeros();

    // compute gradient
    vec gradient = m_system->getWaveFunction()->computeGradient(particles, i);

    force = 2*gradient;

    return force;
}

vec EllipticalHarmonicOscillator::computeGradientPsi(std::vector<Particle*> particles){
    // compute the gradient of wavefunction psi

    int nDim, nPart;
    double psi;
    nDim = m_system->getNumberOfDimensions();
    nPart = m_system->getNumberOfParticles();

    psi = m_system->getWaveFunctionValue();

    vec gradient(nDim);

    gradient = 0.0;
    for (int i=0; i<nPart; i++){
        gradient += m_system->getWaveFunction()->computeGradient(particles, i) * psi;
    }

    return gradient;
}
