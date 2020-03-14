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

double SimpleGaussian::evaluate(std::vector<class Particle*> particles, double alpha) {
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i].getPosition()
     * function.
     *
     * For the actual expression, use exp(-alpha * r^2), with alpha being the
     * (only) variational parameter.
     */
//    vec psi(particles.size()); psi.zeros();
    double psi=1;
    double r;
    int dim = particles[0]->getPosition().size();

    for(int i=0; i<particles.size(); i++){
        for(int j=0; j<dim; j++){
            r = particles[i]->getPosition()[j];     //getPosition()[dimension]
            psi = psi*exp(-alpha*r*r);
        }

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
       h = 1e-7;
       int dim = m_system->getNumberOfDimensions();
       int nPart = m_system->getNumberOfParticles();

       double alpha = m_system->getWaveFunction()->getParameters()[0]; //function getParameters is in class Wavefunction. Wavefunction can be accessed through m_system, as it is defined in m_system

//       cout << "Before loop: " << "\n";
//       cout << particles[0]->getPosition()[0] << endl;

       for (int i = 0; i<nPart; i++){
           for (int j = 0; j<dim; j++){
               // New position
               particles[i]->adjustPosition(h, j);
               wfnew = evaluate(particles,alpha);
//               cout << "New: " << particles[i]->getPosition()[j] << endl;

               // Old position
               particles[i]->adjustPosition(-2*h, j);
               wfold = evaluate(particles, alpha);
//               cout << "Old: " << particles[i]->getPosition()[j] << endl;

               // Current position
               particles[i]->adjustPosition(h, j);
               wfcurrent = evaluate(particles, alpha);
//               cout << "Current: " << particles[i]->getPosition()[j] << endl;

               deriv += (wfnew -2*wfcurrent + wfold)/(h*h);
           }
       }


//       cout << "After loop: " << "\n";
//       cout << particles[0]->getPosition()[0] << endl;

       return deriv;
}
