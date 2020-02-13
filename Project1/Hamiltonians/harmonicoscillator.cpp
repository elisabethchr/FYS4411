#include "harmonicoscillator.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

HarmonicOscillator::HarmonicOscillator(System* system, double omega) :
        Hamiltonian(system) {
    assert(omega > 0);
    m_omega  = omega;
}

double HarmonicOscillator::computeLocalEnergy(std::vector<Particle*> particles) {

    /* Here, you need to compute the kinetic and potential energies. Note that
     * when using numerical differentiation, the computation of the kinetic
     * energy becomes the same for all Hamiltonians, and thus the code for
     * doing this should be moved up to the super-class, Hamiltonian.
     *
     * You may access the wave function currently used through the
     * getWaveFunction method in the m_system object in the super-class, i.e.
     * m_system->getWaveFunction()...
     */

    double h = m_system->m_stepLength;

    double psi_current = m_system->getWaveFunction()->evaluate(particles, 0.4);//how to get alpha?
    std::vector<particle*> particles_new = particles;

    for (int i = 0; i<particles_new.size();i++){
         particles_new->adjust_position(h,0);
    }

    double psi_new = m_system->getWaveFunction()->evaluate(particles_new,0.4);


    particles_new->


    int dim = m_system->getNumberOfDimensions();
    int nPart= m_system->getNumberOfParticles();




    double El = dim*0.5*m_omega*nPart;
    double potentialEnergy = 0;
    double kineticEnergy   = 0;
    return kineticEnergy + potentialEnergy;
}

