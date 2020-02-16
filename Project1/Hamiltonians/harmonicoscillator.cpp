#include "harmonicoscillator.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;
using namespace arma;

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

    int dim = m_system->getNumberOfDimensions();
    int nPart = m_system->getNumberOfParticles();

    // Reset energy for each individual particle (but with renewed wavefunctions)
    double El = dim*0.5*m_omega*nPart;      //analytic expression
    double potentialEnergy = 0;
    double kineticEnergy   = 0;

    // Calculate kinetic energy:
    kineticEnergy -= m_system->getWaveFunction()->computeDoubleDerivative(particles);
    cout << "kinetic energy: " << kineticEnergy << endl;

    return kineticEnergy + potentialEnergy;
}

