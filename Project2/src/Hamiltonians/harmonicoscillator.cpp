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

    double kineticEnergy = 0;
    double potentialEnergy = 0;
    double interactionEnergy = 0;
    int dim = m_system->getNumberDimensions();
    int nPart= m_system->getNumberParticles();
    int M = m_system->getNumberVisibleNodes();



//    arma::vec X = m_system->getWaveFunction()->get_X(); //Gets visible nodes
//    double psi = m_system->getWaveFunction()->evaluate(X); //Gets the wavefunction evaluated at the current pos.

    kineticEnergy = m_system->getWaveFunction()->computeDoubleDerivative_analytic();

//    int particleCounter = 0;
    for(int i = 0; i<M; i++){
        potentialEnergy += m_omega*m_omega*X[i]*X[i];
    }


//Uncomment to include interaction part:

//    for (int p=0; p<nPart-1; p++){
//        for(int q=p+1; q<nPart; q++){
//            double distance = m_system->getWaveFunction()->getDistance(p,q);
//            interactionEnergy+= 1/distance;
//        }
//    }

    return kineticEnergy + potentialEnergy + interactionEnergy;
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
