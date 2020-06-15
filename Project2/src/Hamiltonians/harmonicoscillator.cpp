#include "harmonicoscillator.h"
#include <cassert>
#include <iostream>
#include <vector>
#include <armadillo>
#include "../system.h"
#include "../WaveFunctions/wavefunction.h"

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
    int nPart = m_system->getNumberParticles();
    int M = m_system->getNumberVisibleNodes();

    //    arma::vec X = m_system->getWaveFunction()->get_X(); //Gets visible nodes
    //    double psi = m_system->getWaveFunction()->evaluate(X); //Gets the wavefunction evaluated at the current pos.

    kineticEnergy = m_system->getWaveFunction()->computeDoubleDerivative_analytic();

    //    int particleCounter = 0;
    for(int i = 0; i<M; i++){
        potentialEnergy += 0.5*m_omega*m_omega*X[i]*X[i];
    }

    // add interaction term if stated so
    if (m_interaction == true){
        for (int p=0; p<nPart-1; p++){
            for(int q=p+1; q<nPart; q++){
                double distance = m_system->getWaveFunction()->getDistance(p,q);
                interactionEnergy+= 1.0/distance;
            }
        }
    }

    arma::vec pos = m_system->getWaveFunction()->get_X();
//    cout << "\nkineticEnergy + potentialEnergy = " << kineticEnergy + potentialEnergy << endl;
//    cout << "interaction energy = " << interactionEnergy << endl;
//    cout << "interaction energy test = " << interaction(pos, M, dim) << endl;

//    cout << "Local energy = " << kineticEnergy + potentialEnergy + interactionEnergy << endl;
    return kineticEnergy + potentialEnergy + interactionEnergy;
}

/*
double HarmonicOscillator::interaction(arma::vec x, int nx, int dim) {
    double interactionTerm = 0;
    double rDistance;
    // Loop over each particle
    for (int r=0; r<nx-dim; r+=dim) {
        // Loop over each particle s that particle r hasn't been paired with
        for (int s=(r+dim); s<nx; s+=dim) {
            rDistance = 0;
            // Loop over each dimension
            for (int i=0; i<dim; i++) {
                rDistance += (x(r+i) - x(s+i))*(x(r+i) - x(s+i));
            }
            interactionTerm += 1.0/sqrt(rDistance);
        }
    }
    return interactionTerm;
}
*/


/* Compute the gradient of the local energy wrt. the RBM parameters, i.e. (1/psi)*d(psi)/d(alpha_i) */
arma::vec HarmonicOscillator::computeLocalEnergyGradient(){
    // get necessary variables
    arma::vec m_a = m_system->getWaveFunction()->get_a();
    arma::vec m_b = m_system->getWaveFunction()->get_b();
    arma::mat m_w = m_system->getWaveFunction()->get_w();
    arma::vec m_x = m_system->getWaveFunction()->get_X();
    m_sigma = m_system->getWaveFunction()->getSigma();
    m_sigma2 = m_sigma*m_sigma;

    m_nv = m_system->getNumberVisibleNodes();
    m_nh = m_system->getNumberHiddenNodes();

    arma::vec O = m_b + (m_x.t()*m_w).t()*(1/((double) m_sigma*m_sigma));
    arma::vec dPsi; dPsi.zeros(m_nv + m_nh + m_nv*m_nh);

    // compute d(psi)/d(a)*1/psi
    for (int i=0; i<m_nv; i++){
        dPsi[i] = (m_x[i] - m_a[i])/(m_sigma2);
    }

    // compute d(psi)/d(b)*1/psi
    for (int i=m_nv; i<m_nv+m_nh; i++){
        dPsi[i] = 1.0/(exp(-O[i - m_nv]) + 1);
    }

    // compute d(psi)/d(w)*1/psi
    int i = m_nv + m_nh;
    for (int j=0; j<m_nv; j++){
        for (int k=0; k<m_nh; k++){
            dPsi[i] = m_x[j]/(exp(-O[k]) + 1)*(1.0/m_sigma2);
            i++;
        }
    }

    return dPsi;
}


/* Compute the drift force experienced by particles */
vec HarmonicOscillator::computeQuantumForce(arma::vec X, int i){

    arma::vec a = zeros(1);

    return a;
}
