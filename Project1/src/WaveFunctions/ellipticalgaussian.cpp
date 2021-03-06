#include "ellipticalgaussian.h"
#include <cmath>
#include <cassert>
#include <armadillo>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

using namespace std;
using namespace arma;

// Obtain parameters required for the evalutaion of the elliptical wavefunction
EllipticalGaussian::EllipticalGaussian(System* system, std::vector<double> alpha, std::vector<double> beta, double hard_core_diameter, bool Jastrow, double gamma) :
    WaveFunction(system) {
    assert(alpha.size() >= 0);

    m_parameters = alpha;
    m_numberOfParameters = alpha.size();

    m_beta = beta;
    m_numberOfParametersBeta = beta.size();

    m_hard_core_diameter = hard_core_diameter;
    m_Jastrow = Jastrow;
    m_gamma = gamma;
}

/*
 * The following expressions are used to evaluate the wavefunction,
 * when considering an elliptical potential along with correlation
 * between the particles
*/

double EllipticalGaussian::evaluate(std::vector<class Particle*> particles) {

    // Evaluate the full expression for the elliptical wavefunction using the single-particle wavefunction,
    // the correlation wavefunction, and the distance between particles.

    double r;
    double psi = 1.0;
    int nDim = m_system->getNumberOfDimensions();
    int nPart = m_system->getNumberOfParticles();

    std::vector<double> alphas = m_parameters;

    for (int i=0; i<nPart; i++){

        psi *= SingleParticleFunction(particles, i);

        for(int j=i+1; j<nPart; j++){
            r = getDistance(particles, i, j);
            psi *= correlationWaveFunction(particles, r);
        }
    }

    return psi;
}


double EllipticalGaussian::SingleParticleFunction(std::vector<class Particle *> particles, int particle) {
    // Compute single-particle function g(alpha, beta, r)

    int i = m_system->getAlphaIndex();
    int nDim = m_system->getNumberOfDimensions();
    double r = 0.0;
    double alpha = m_parameters[i];
    double beta = m_beta[0];

    std::vector<double> pos = particles[particle]->getPosition();

    for(int j=0; j<nDim; j++){

        // perturb the z-direction by beta
        if(j==2){ r += beta*pos[j]*pos[j]; }
        else { r += pos[j]*pos[j]; }
    }

    return exp(-alpha*r);
}


double EllipticalGaussian::correlationWaveFunction(std::vector<class Particle *> particles, double distance) {
    // Compute correlation wave function f(a, r), i.e. Jastrow factor

    double a = m_hard_core_diameter;
    double f = 1;

    if (distance <= a){ f = 0; }

    else {f = 1 - a / ((double) distance); }

    if(!m_Jastrow){f=1;}

    return f;
}


double EllipticalGaussian::getDistance(std::vector<class Particle *> particles, int particle_i, int particle_j) {
    // Obtain the distance between two particles

    int nDim = m_system->getNumberOfDimensions();
    double distance = 0.0;
    double diff = 0.0;

    std::vector<double> pos_i = particles[particle_i]->getPosition();
    std::vector<double> pos_j = particles[particle_j]->getPosition();

    for (int j=0; j<nDim; j++){
        diff = abs(pos_i[j] - pos_j[j]);
        distance += diff*diff;
    }

    return sqrt(distance);
}



/*
 * The following functions are used to compute the analytical expression for the laplacian,
 * when including the Jastrow factor
*/

double EllipticalGaussian::computeLaplacian(std::vector<class Particle *> particles){
    int m_numberOfParticles = m_system->getNumberOfParticles();
    int m_numberOfDimensions = m_system->getNumberOfDimensions();

    arma::vec gradient_SPF(m_numberOfDimensions);
    arma::vec gradient_correlation(m_numberOfDimensions);

    double laplacian = 0.0;

    double term1, term2, term3, term4;
    term1 = 0.0; term2 = 0.0; term3 = 0.0; term4 = 0.0;

    for (int k=0; k<m_numberOfParticles; k++){
        term1 += computeLaplacian_SPF(k);
        term4 += computeLaplacian_correlation(particles, k);

        gradient_SPF = computeGradient_SPF(k);
        gradient_correlation = computeGradient_correlation(particles, k);
        for (int dim=0; dim<m_numberOfDimensions; dim++){
            term2 += 2*gradient_SPF[dim] * computeLaplacian_correlation(particles, k);
            term3 += gradient_correlation[dim]*gradient_correlation[dim];
        }
    }
    laplacian += term1 + term2 + term3 + term4;

    return laplacian;
}


arma::vec EllipticalGaussian::computeGradient(std::vector<class Particle *> particles, int i) {
    // compute the gradient

    int nDim, nPart;
    double term1, term2, rij, sum1, sum2;
    nDim = m_system->getNumberOfDimensions();
    nPart = m_system->getNumberOfParticles();

    arma::vec gradient(nDim);
    arma::vec gradient_SPF(nDim);
    arma::vec gradient_correlation(nDim);

    term1 = 0.0; term2 = 0.0; sum1 = 0.0; sum2 = 0.0; gradient = 0.0;

    gradient_SPF = computeGradient_SPF(i);
    gradient_correlation = computeGradient_correlation(particles, i);

    for (int dim=0; dim<nDim; dim++){
        gradient[dim] += gradient_SPF[dim] + gradient_correlation[dim];
    }
    return gradient;
}


double EllipticalGaussian::computeDerivativePsi_alpha(std::vector<class Particle *> particles){
    // compute the derivative of psi wrt. variational parameter alpha

    int nPart, nDim;
    double m_alpha, m_beta, deriv, r, rij, psi;
    std::vector<double> pos;

    nDim = m_system->getNumberOfDimensions();
    nPart = m_system->getNumberOfParticles();
    m_alpha = getParameters()[0];
    m_beta = getParametersBeta()[0];
//    psi = evaluate(particles);
    psi = m_system->getWaveFunctionValue();

    for (int i=0; i<nPart; i++){
        pos = particles[i]->getPosition();
        for (int dim=0; dim<nDim; dim++){
            r = pos[dim]*pos[dim];
            if (dim==2){ deriv += m_beta*r; }
            else { deriv += r; }
        }
    }

    deriv *= psi;

    return deriv;
}


double EllipticalGaussian::computeLaplacian_SPF(int particle_i){
    // compute the laplacian of the single-particle wavefunction

    double r, laplacian;
    int nDim = m_system->getNumberOfDimensions();

    int i = m_system->getAlphaIndex();
    double alpha = getParameters()[i];
    double beta = getParametersBeta()[0];

    std::vector<double> pos = m_system->getParticles()[particle_i]->getPosition();

    laplacian = 0.0;

    for (int j=0; j<nDim; j++){
        r = pos[j];
        // perturb the z-direction by beta
        if (j==2){ r *= beta; }

        laplacian += -2*alpha*r*r;
    }

    laplacian *= -2*alpha;

    if (nDim ==3) {laplacian += -2*alpha*(nDim - 1 + beta); }
    else { laplacian += -2*alpha*(nDim); }

    return laplacian;
}


vec EllipticalGaussian::computeGradient_SPF(int particle_i){
    // compute the gradient of the single-particle wavefunction

    double r;
    int nDim = m_system->getNumberOfDimensions();

    int i = m_system->getAlphaIndex();
    double alpha = getParameters()[i];
    double beta = getParametersBeta()[0];

    arma::vec gradient(nDim);
    std::vector<double> pos = m_system->getParticles()[particle_i]->getPosition();

    for(int j=0; j<nDim; j++){
        r = pos[j];
        if(nDim==3){ gradient[j] = -2*alpha*beta*r; }
        else { gradient[j] = -2*alpha*r; }
    }

    return gradient;
}


double EllipticalGaussian::computeLaplacian_correlation(std::vector<class Particle *> particles, int particle_i){
    // compute the laplacian of the correlation wavefunction

    double a, rij, du, ddu, laplacian;
    int nPart = m_system->getNumberOfParticles();

    a = m_hard_core_diameter;

    // initialise position vector for particle i
    std::vector<double> pos_i;
    pos_i = m_system->getParticles()[particle_i]->getPosition();

    laplacian = 0.0;
    for (int particle_j=0; particle_i<nPart; particle_i++){
        if (particle_j != particle_i){
            rij = getDistance(particles, particle_i, particle_j);
            // ensure boundary conditions for hard-sphere diameter
            if (rij <= a){ laplacian = 0; }
            else{
                du = computeDerivative_correlation(particle_i, particle_j);
                ddu = computeDoubleDerivative_correlation(particle_i, particle_j);

                laplacian += ddu + 2*du/rij;
            }
        }
    }

    return laplacian;
}


arma::vec EllipticalGaussian::computeGradient_correlation(std::vector<class Particle *> particles, int particle_i){
    // compute the gradient of the correlation wavefunction

    double a = m_hard_core_diameter;
    double rij, du;
    int nPart = m_system->getNumberOfParticles();
    int nDim = m_system->getNumberOfDimensions();

    std::vector<double> pos_i;
    std::vector<double> pos_j;

    arma::vec gradient(nDim);

    pos_i = m_system->getParticles()[particle_i]->getPosition();

    gradient = 0.0;

    for (int particle_j=0; particle_j<nPart; particle_j++){

        if (particle_i != particle_j){
            rij = getDistance(particles, particle_i, particle_j);
            pos_j = m_system->getParticles()[particle_j]->getPosition();

            du = computeDerivative_correlation(particle_i, particle_j);

            for (int dim=0; dim<nDim; dim++){
                // ensure boundary conditions for hard-sphere diameter
                if (rij <= a){ gradient[dim] = 0; }
                else{
                    gradient[dim] += (pos_i[dim] - pos_j[dim]) /rij * du;
                }
            }
        }
    }
    return gradient;
}


double EllipticalGaussian::computeDerivative_correlation(int particle_i, int particle_j){
    // compute the derivative of the correlation functoin wrt. distance between particles (rij)

    double a = m_hard_core_diameter;
    int nDim = m_system->getNumberOfDimensions();

    std::vector<class Particle* > particles = m_system->getParticles();

    double deriv = 0.0;
    double rij = getDistance(particles, particle_i, particle_j);

    deriv -= a / (rij*(rij - a));

    return deriv;
}


double EllipticalGaussian::computeDoubleDerivative_correlation(int particle_i, int particle_j){
    // compute the double derivative of the correlation function wrt. distance between particles (rij)

    double a = m_hard_core_diameter;
    double deriv = 0.0;

    std::vector<class Particle* > particles = m_system->getParticles();

    double rij = getDistance(particles, particle_i, particle_j);

    deriv += a*(a-2*rij) / ((double) rij*rij*(a - rij)*(a - rij));

    return deriv;
}


double EllipticalGaussian::computeDoubleDerivative_analytic(std::vector<class Particle*> particles) {
     // Compute the analytic laplacian used in finding the expression for the analytic local energy

    return computeLaplacian(particles);
}


double EllipticalGaussian::computeDoubleDerivative_numeric(std::vector<class Particle*> particles) {
    // Compute the numeric laplacian used in finding the expresion for the numeric local energy

    double h, wfplus, wfminus, wfcurrent;
    double deriv;
    int dim = m_system->getNumberOfDimensions();
    int nPart = m_system->getNumberOfParticles();

    // steplength
    h = 1e-4;

    // compute wavefunction for current positions
    wfcurrent = m_system->getWaveFunctionValue();

    // perform the numeric double derivative
    for (int i = 0; i<nPart; i++){
        for (int j = 0; j<dim; j++){
            // Change in positive direction
            particles[i]->adjustPosition(h, j);
            wfplus = evaluate(particles);
            //               cout << "New: " << particles[i]->getPosition()[j] << endl;

            // Change in negative direction
            particles[i]->adjustPosition(-2*h, j);
            wfminus = evaluate(particles);
            //               cout << "Old: " << particles[i]->getPosition()[j] << endl;

            // Current position
            particles[i]->adjustPosition(h, j);

            deriv += (wfplus -2*wfcurrent + wfminus)/(h*h);
        }
    }
    return deriv;
}


