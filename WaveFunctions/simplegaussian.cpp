#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include <armadillo>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

using namespace std;
using namespace arma;

SimpleGaussian::SimpleGaussian(System* system, std::vector<double> alpha) :
    WaveFunction(system) {
    assert(alpha.size() >= 0);

    //    m_numberOfParameters = 1;
    //    m_parameters.reserve(1);
    //    m_parameters.push_back(alpha);
    m_parameters = alpha;
    m_numberOfParameters = alpha.size();
}

/*
SimpleGaussian::SimpleGaussian(System* system, double alpha) :
    WaveFunction(system) {
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
}
*/

double SimpleGaussian::evaluate(std::vector<class Particle*> particles) {
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
    std::vector<double> alphas = m_parameters;

    //    for(const auto &i:alphas){cout << "alphas = " << i << endl;}

    int i = m_system->getAlphaIndex();
    double alpha = m_parameters[i];
    //    if(m_system->getMetropolisStep() == 0){ cout << "alpha used: " << alpha << endl; }

    for(int i=0; i<particles.size(); i++){
        for(int j=0; j<dim; j++){
            r = particles[i]->getPosition()[j];     //getPosition()[dimension]
            psi = psi*exp(-alpha*r*r);
        }

    }
    return psi;
}

arma::vec SimpleGaussian::computeGradient(std::vector<class Particle *> particles, int i){
    double h, wfnew, wfold;
    h = 1e-4;
    int dim = m_system->getNumberOfDimensions();
//    int nPart = m_system->getNumberOfParticles();

    //    double alpha = m_system->getWaveFunction()->getParameters()[0]; //function getParameters is in class Wavefunction. Wavefunction can be accessed through m_system, as it is defined in m_system

    arma::vec gradient(dim); gradient.zeros();

    wfold = m_system->getWaveFunction()->evaluate(particles);

    for (int j = 0; j<dim; j++){
        // New position
        particles[i]->adjustPosition(h, j);
        wfnew = evaluate(particles);

        // Current position
        particles[i]->adjustPosition(-h, j);

        gradient[j] = (wfnew - wfold)/(h);
    }

    return gradient;
}

double SimpleGaussian::computeDoubleDerivative_numeric(std::vector<class Particle*> particles) {
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
    //h = 1e-7;
    h = 1e-4;
    int dim = m_system->getNumberOfDimensions();
    int nPart = m_system->getNumberOfParticles();

    int i = m_system->getAlphaIndex();
    double alpha = m_parameters[i];

    wfcurrent = m_system->getWaveFunctionValue();


    //    double alpha = m_system->getWaveFunction()->getParameters()[0]; //function getParameters is in class Wavefunction. Wavefunction can be accessed through m_system, as it is defined in m_system

    //       cout << "Before loop: " << "\n";
    //       cout << particles[0]->getPosition()[0] << endl;

    for (int i = 0; i<nPart; i++){
        for (int j = 0; j<dim; j++){
            // New position
            particles[i]->adjustPosition(h, j);
            wfnew = evaluate(particles);
            //               cout << "New: " << particles[i]->getPosition()[j] << endl;

            // Old position
            particles[i]->adjustPosition(-2*h, j);
            wfold = evaluate(particles);
            //               cout << "Old: " << particles[i]->getPosition()[j] << endl;

            // Current position
            particles[i]->adjustPosition(h, j);
            //wfcurrent = evaluate(particles);
            //            cout << "wfcurrent_getWFValue = " << wfcurrent_test << endl;
            //            cout << "wfcurrent = " << wfcurrent << endl;
            //               cout << "Current: " << particles[i]->getPosition()[j] << endl;

            deriv += (wfnew -2*wfcurrent + wfold)/(h*h);
        }
    }
    return deriv;
}

double SimpleGaussian::computeDoubleDerivative_analytic(std::vector<class Particle*> particles) {
    /* All wave functions need to implement this function, so you need to
        * find the double derivative analytically. Note that by double derivative,
        * we actually mean the sum of the Laplacians with respect to the
        * coordinates of each particle.
        *
        * This quantity is needed to compute the (local) energy (consider the
        * Schrödinger equation to see how the two are related).
        */
    double pos, deriv;
    int dim = m_system->getNumberOfDimensions();
    int nPart = m_system->getNumberOfParticles();

    int i = m_system->getAlphaIndex();
    double alpha = m_parameters[i];

    double wf = m_system->getWaveFunctionValue();//evaluate(particles);

    //    double alpha = m_system->getWaveFunction()->getParameters()[0]; //function getParameters is in class Wavefunction. Wavefunction can be accessed through m_system, as it is defined in m_system

    //       cout << "Before loop: " << "\n";
    //       cout << particles[0]->getPosition()[0] << endl;

    for (int i = 0; i<nPart; i++){
        for (int j = 0; j<dim; j++){
            pos = particles[i]->getPosition()[j];
            deriv -= 2*alpha*(2*alpha*pos*pos)*wf;
        }
        deriv += 2*alpha*dim*wf;
    }
    return deriv;
}
