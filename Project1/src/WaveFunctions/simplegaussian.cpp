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

    m_parameters = alpha;
    m_numberOfParameters = alpha.size();
}

double SimpleGaussian::evaluate(std::vector<class Particle*> particles) {

    double psi=1;
    double r;
    int dim = particles[0]->getPosition().size();
    std::vector<double> alphas = m_parameters;

    int i = m_system->getAlphaIndex();
    double alpha = m_parameters[i];

    for(int i=0; i<particles.size(); i++){
        for(int j=0; j<dim; j++){
            r = particles[i]->getPosition()[j];
            psi = psi*exp(-alpha*r*r);
        }

    }
    return psi;
}

arma::vec SimpleGaussian::computeGradient(std::vector<class Particle *> particles, int i){
    double h, wfnew, wfold;
    h = 1e-4;
    int dim = m_system->getNumberOfDimensions();

    arma::vec gradient(dim); gradient.zeros();

//    wfold = m_system->getWaveFunction()->evaluate(particles);
    wfold = m_system->getWaveFunctionValue();

    for (int j = 0; j<dim; j++){
        // New position
        particles[i]->adjustPosition(h, j);
        wfnew = evaluate(particles);

        // Current position
        particles[i]->adjustPosition(-h, j);

        gradient[j] = (wfnew - wfold)/h;
    }

    return gradient;
}

double SimpleGaussian::computeDoubleDerivative_numeric(std::vector<class Particle*> particles) {

    double h, wfnew, wfold, wfcurrent;
    double deriv;
    //h = 1e-7;
    h = 1e-4;
    int dim = m_system->getNumberOfDimensions();
    int nPart = m_system->getNumberOfParticles();

    wfcurrent = m_system->getWaveFunctionValue();

    for (int i = 0; i<nPart; i++){
        for (int j = 0; j<dim; j++){
            // New position
            particles[i]->adjustPosition(h, j);
            wfnew = evaluate(particles);

            // Old position
            particles[i]->adjustPosition(-2*h, j);
            wfold = evaluate(particles);

            // Current position
            particles[i]->adjustPosition(h, j);

            deriv += (wfnew -2*wfcurrent + wfold)/(h*h);
        }
    }
    return deriv;
}

double SimpleGaussian::computeDoubleDerivative_analytic(std::vector<class Particle*> particles) {

    double pos, deriv;
    int dim = m_system->getNumberOfDimensions();
    int nPart = m_system->getNumberOfParticles();

    int i = m_system->getAlphaIndex();
    double alpha = m_parameters[i];

    double wf = m_system->getWaveFunctionValue();//evaluate(particles);

    for (int i = 0; i<nPart; i++){
        for (int j = 0; j<dim; j++){
            pos = particles[i]->getPosition()[j];
            deriv -= 2*alpha*(2*alpha*pos*pos)*wf;
        }
        deriv += 2*alpha*dim*wf;
    }
    return deriv;
}
