#include <cmath>
#include <cassert>
#include <armadillo>
#include <random>
#include "neuralquantumstate.h"
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

using namespace std;
using namespace arma;

NeuralQuantumState::NeuralQuantumState(System* system, int n_hidden, int n_visible, int part, int dim, double sigma) :
    WaveFunction(system) {
    assert(n_hidden > 0 && n_visible > 0);
    m_nh = n_hidden;
    m_nv = n_visible;
    m_part = part;
    m_dim = dim;
    m_sigma = sigma;

    std::random_device rd;
    m_randomEngine = std::mt19937_64(rd());

    m_system->setNumberHiddenNodes(n_hidden);
    m_system->setNumberVisibleNodes(n_visible);
    m_system->setNumberParticles(part);
    m_system->setNumberDimensions(dim);

    setupInitialState();
}

/* Set up initial position, hidden and visible biases, and interaction between hidden and visible biases */
void NeuralQuantumState::setupInitialState(){
    m_x.zeros(m_nv);
    m_a.zeros(m_nv);
    m_b.zeros(m_nh);
    m_w.zeros(m_nv, m_nh);

    std::uniform_real_distribution<double> uniform_weights(-1.0,1.0);
    std::uniform_real_distribution<double> uniform_position(-0.5,0.5);

    for (int i=0; i<m_nv; i++){
        m_x[i] = uniform_position(m_randomEngine);
        m_a[i] = uniform_weights(m_randomEngine);
        for (int j=0; j<m_nh; j++){
            m_w(i, j) = uniform_weights(m_randomEngine);
        }
    }

    for (int j=0; j<m_nh; j++){
        m_b[j] = uniform_weights(m_randomEngine);
    }
}


/* Evaluate wavefunction at a certain position */
double NeuralQuantumState::evaluate(arma::vec position){
    m_term1 = 0.0;
    for (int i=0; i<m_nv; i++){
        m_term1 -= (m_x[i] - m_a[i])*(m_x[i] - m_a[i]);
    }
    m_term1 = m_term1/(2*m_sigma*m_sigma);

    m_term2 = 1.0;
    O = m_b + ((m_x.t()*m_w).t()) / ((double) m_sigma*m_sigma);
    for (int j=0; j<m_nh; j++){
        m_term2 *= (1 + exp(O[j]));
    }

    return exp(m_term1)*m_term2;
}

double NeuralQuantumState::sigmoid(double x){
    return (1/(1+std::exp(-x));
}

double NeuralQuantumState::v(int j, arma::vec m_x, arma::mat m_w){

    return m_b[j] + m_x*m_w.col(j);
}

arma::mat NeuralQuantumState::hadamardProd(arma::mat w){
    int M = m_system->getNumberVisibleNodes();
    int N = m_system->getNumberHiddenNodes();

    arma::mat w_Hadard;
    w_Hadard.zeros(M,N);
    for(int i; i<M; i++){
        for(int j; j<N; j++){
            w_Hadamard(i,j) = w(i,j)*w(i,j);
        }
    }
    return w_Hadarad;
}

double NeuralQuantumState::getDistance(int p, int q){
    int dim = m_system->getNumberDimensions();

    double dist = 0;



    for(int d = 0; d<dim; d++){
        dist+= (X[dim*p+d]-X[dim*q+d])*(X[dim*p+d]-X[dim*q+d]); //Assumes particle index starts at 0
    }

    return sqrt(dist);

}

double NeuralQuantumState::computeDoubleDerivative_analytic(){

    int M = m_system->getNumberVisibleNodes();
    int N = m_system->getNumberHiddenNodes();

    arma::vec S;
    arma::vec S_tilde;
    arma::vec one_vector;
    S.zeros(N);
    S_tilde.zeros(N;
    one_vector.ones(M);

    for (int j = 0; j<N; j++){
        S[j] = sigmoid(v(j,m_x, m_w));
        S_tilde[j] = sigmoid(v(j,m_x, m_w))*sigmoid(-v(j,m_x, m_w));
    }


    double E_K= -1/(2*pow(m_sigma,4))*((m_x-m_a)*(m_x-m_a).t() - 2*S*(m_w.t()*(m_x-m_a).t())+(S*m_w.t())*(m_w*S.t())+one_vector*(hadamardProd(m_w)*S_tilde.t())) - M*(2*pow(sigma,2));

    return E_k;
}


