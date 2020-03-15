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

    /* Here, you need to compute the kinetic and potential energies. Note that
     * when using numerical differentiation, the computation of the kinetic
     * energy becomes the same for all Hamiltonians, and thus the code for
     * doing this should be moved up to the super-class, Hamiltonian.
     *
     * You may access the wave function currently used through the
     * getWaveFunction method in the m_system object in the super-class, i.e.
     * m_system->getWaveFunction()...
     */

    double potentialEnergy = 0;
    double kineticEnergy   = 0;
    int dim = m_system->getNumberOfDimensions();
    int nPart= m_system->getNumberOfParticles();
    //    std::vector<double> alpha = m_system->getWaveFunction()->getParameters();

    for (int i = 0; i<nPart;i++){
        for(int j = 0; j<dim; j++){
            double coordinate = particles.at(i)->getPosition().at(j);
            potentialEnergy+= 0.5*m_omega*m_omega*coordinate*coordinate;
        }
    }

    // possible to setWavefunction in MetropolisStep and avoid extra evaluate? (setWaveFunction(psi), m_psi = getWavefunction->getValue();)
    //    double psi = m_system->getWaveFunction()->evaluate(particles);

    double psi = m_system->getWaveFunctionValue();

    if(type == true){ kineticEnergy -= (1/psi)*0.5*m_system->getWaveFunction()->computeDoubleDerivative_numeric(particles); }
    else if(type == false){ kineticEnergy = (1/psi)*0.5*m_system->getWaveFunction()->computeDoubleDerivative_analytic(particles); }

    double El = potentialEnergy + kineticEnergy;
    //    cout <<"Local energy:  " << El << endl;
    return El;
}

vec HarmonicOscillator::computeQuantumForce(std::vector<Particle *> particles, int i){
//    double h = 1e-7;
    int nDim = m_system->getNumberOfDimensions();
//    int nPart = m_system->getNumberOfParticles();
    vec force(nDim); force.zeros();
    vec gradient = m_system->getWaveFunction()->computeGradient(particles, i);

    force = 2*gradient;

//    cout << "computeQuantumForce okay" << endl;

    return force;
}
/*
//    for (int i=0; i<nPart; i++){
        for (int j=0; j<dim; j++){
            // position in positive direction
            particles[i]->adjustPosition(h, j);
            wfNew = m_system->getWaveFunction()->evaluate(particles);
            // position in negative direction
            particles[i]->adjustPosition(-h, j);
//            wfminus = m_system->getWaveFunction()->evaluate(particles);
            // calculate derivative
            deriv(i, j) = (wfNew - wfOld)/h;
//        }
    }
*/
//    return deriv;
//}










