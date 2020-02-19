#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include <armadillo>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

using namespace std;
using namespace arma;

SimpleGaussian::SimpleGaussian(System* system, double alpha) :
    WaveFunction(system) {
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
}

double SimpleGaussian::evaluate(std::vector<class Particle*> particles) {
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i].getPosition()
     * function.
     *
     * For the actual expression, use exp(-alpha * r^2), with alpha being the
     * (only) variational parameter.
     */
    double psi = 1;
    double r;
    double alpha = m_parameters[0];

    for(int i=0; i<particles.size(); i++){
        r = particles[i]->getPosition()[0];     //getPosition()[dimension]
        psi *= exp(-alpha*r*r);
    }
    return psi;
}

double SimpleGaussian::computeDoubleDerivative(std::vector<class Particle*> particles) {
    /* All wave functions need to implement this function, so you need to
     * find the double derivative analytically. Note that by double derivative,
     * we actually mean the sum of the Laplacians with respect to the
     * coordinates of each particle.
     *
     * This quantity is needed to compute the (local) energy (consider the
     * Schrödinger equation to see how the two are related).
     */
    double h, wfnew, wfold, wfcurrent;
    double deriv;
    h = m_system->getStepLength();
    int dim = m_system->getNumberOfDimensions();
    int nPart = m_system->getNumberOfParticles();

    double alpha = m_system->getWaveFunction()->getParameters()[0]; //function getParameters is in class Wavefunction. Wavefunction can be accessed through m_system, as it is defined in m_system

    cout << "Before loop: " << "\n";
    cout << particles[0]->getPosition()[0] << endl;

    for (int i = 0; i<nPart; i++){
        for (int j = 0; j<dim; j++){
            // New position
            particles[i]->adjustPosition(h, j);
            wfnew = evaluate(particles);
            cout << "New: " << particles[i]->getPosition()[j] << endl;

            // Old position
            particles[i]->adjustPosition(-2*h, j);
            wfold = evaluate(particles);
            cout << "Old: " << particles[i]->getPosition()[j] << endl;

            // Current position
            particles[i]->adjustPosition(h, j);
            wfcurrent = evaluate(particles);
            cout << "Current: " << particles[i]->getPosition()[j] << endl;
        }
    }

    cout << "After loop: " << "\n";
    cout << particles[0]->getPosition()[0] << endl;

    deriv = (wfnew + wfcurrent - 2*wfold)/(h*h);

    return deriv;
}
