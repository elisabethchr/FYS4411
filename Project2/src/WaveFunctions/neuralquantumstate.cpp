#include <cmath>
#include <cassert>
#include <armadillo>
#include <random>
#include "neuralquantumstate.h"
#include "wavefunction.h"
#include "../system.h"

using namespace std;
using namespace arma;

NeuralQuantumState::NeuralQuantumState(System* system, int n_hidden, int n_visible, int part, int dim, double sigma, bool gaussian, double initialization) :
    WaveFunction(system) {
    assert(n_hidden > 0 && n_visible > 0);
    m_nh = n_hidden;
    m_nv = n_visible;
    m_part = part;
    m_dim = dim;
    m_sigma = sigma;
    m_gaussian = gaussian;
    m_initialization = initialization;

    std::random_device seed;
    mt19937_64 gen(seed());
    mt19937_64 m_randomengine = gen;

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

//    std::uniform_real_distribution<double> uniform_weights(-1.0 ,1.0); //original
//    std::uniform_real_distribution<double> uniform_weights(-0.10 ,0.10); // new initialization
//    std::uniform_real_distribution<double> uniform_weights(-1e-7, 1e-7); // new initialization 1e-7
    std::uniform_real_distribution<double> uniform_weights(-m_initialization, m_initialization);
    std::uniform_real_distribution<double> uniform_position(-0.5, 0.5);
    std::normal_distribution<double> normal_weights(0, m_initialization);

    // initialize weights according to either a uniform or gaussian distribution
    if (m_gaussian == false){
        for (int i=0; i<m_nv; i++){
            m_a[i] = uniform_weights(m_randomEngine);
            for (int j=0; j<m_nh; j++){
                m_w(i, j) = uniform_weights(m_randomEngine);
            }
        }

        for (int j=0; j<m_nh; j++){
            m_b[j] = uniform_weights(m_randomEngine);
        }
    }

    else if (m_gaussian == true){
        for (int i=0; i<m_nv; i++){
            m_a[i] = normal_weights(m_randomEngine);
            for (int j=0; j<m_nh; j++){
                m_w(i, j) = normal_weights(m_randomEngine);
            }
        }

        for (int j=0; j<m_nh; j++){
            m_b[j] = normal_weights(m_randomEngine);
        }
    }

    for (int i=0; i<m_nv; i++){
        m_x[i] = uniform_position(m_randomEngine);
    }
}


/* Evaluate wavefunction at a certain position */
double NeuralQuantumState::evaluate(arma::vec position){
    m_term1 = 0.0;
    for (int i=0; i<m_nv; i++){
        m_term1 += (position[i] - m_a[i])*(position[i] - m_a[i]);
    }
    m_term1 = exp(-m_term1/(2*m_sigma*m_sigma));

    m_term2 = 1.0;
    O = m_b + ((position.t()*m_w).t()) / ((double) m_sigma*m_sigma);
    for (int j=0; j<m_nh; j++){
        m_term2 *= (1 + exp(O[j]));
    }

    return m_term1*m_term2;
}



/* Compute the logistic function of x */
double NeuralQuantumState::sigmoid(double x){
    return (1/(1+exp(-x)));
}

/* Compute the exponent of the exponential in the logistic function */
double NeuralQuantumState::v(int j){
    double a = (m_b[j] + m_x.t()*m_w.col(j)).eval()(0,0);
    return (m_b[j] + m_x.t()*m_w.col(j)).eval()(0,0);
}

/* Compute the squared of the interaction matrix w */
arma::mat NeuralQuantumState::hadamardProd(arma::mat w){
    int M = m_system->getNumberVisibleNodes();
    int N = m_system->getNumberHiddenNodes();

    arma::mat w_Hadamard;
    w_Hadamard.zeros(M,N);
    for(int i; i<M; i++){
        for(int j; j<N; j++){
            w_Hadamard(i,j) = w(i,j)*w(i,j);
        }
    }
    return w_Hadamard;
}

/* Get distance between particles */
double NeuralQuantumState::getDistance(int p, int q){
    int dim = m_system->getNumberDimensions();
    double dist = 0;

    for(int d = 0; d<dim; d++){
        dist+= (m_x[dim*p+d]-m_x[dim*q+d])*(m_x[dim*p+d]-m_x[dim*q+d]); // assumes particle index starts at 0
    }

    return sqrt(dist);
}

/* Compute the derivative of psi wrt. position coordinates (d(ln psi)/d(X_i)) */
double NeuralQuantumState::computeDerivative_analytic(arma::vec position, int coor){

    arma::vec O = m_b + (m_x.t()*m_w).t()/(m_sigma*m_sigma);

    double deriv=0.0;
    double sum1 = 0.0;
    for(int j=0; j<m_nh; j++){
        sum1 += m_w(coor, j)/(m_sigma*m_sigma*(1.0+exp(-O[j])));
    }
    deriv = -(position[coor]-m_a[coor])/(m_sigma*m_sigma) + sum1/(m_sigma*m_sigma);

    return deriv;
}

/* Compute the double derivative of psi wrt. position coordinates */
double NeuralQuantumState::computeDoubleDerivative_analytic(){

    int M = m_system->getNumberVisibleNodes();
    int N = m_system->getNumberHiddenNodes();

    arma::vec S; S.zeros(N);
    arma::vec S_tilde; S_tilde.zeros(N);
    arma::vec one_vector; one_vector.ones(M);

    for (int j = 0; j<N; j++){
        S[j] = sigmoid(v(j));
        S_tilde[j] = sigmoid(v(j))*sigmoid(-v(j));
    }

    arma::vec O = m_b + (m_x.t()*m_w).t()/(m_sigma*m_sigma);

    double Ek = 0.0;
    double term1 = 0.0;
    double term2 = 0.0;
    double term3 = 0.0;
    double term4 = 0.0;
    for (int i=0; i<m_nv; i++){
        double sum1 = 0.0;
        double sum2 = 0.0;
        double sum3 = 0.0;
        for (int j=0; j<m_nh; j++){
            sum1 += m_w(i, j)*S[j];
            sum2 += m_w(i, j)*S[j];
            sum3 += m_w(i, j)*m_w(i, j)*S_tilde[j];
        }
        term1 += (m_x[i] - m_a[i])*(m_x[i] - m_a[i]);
        term2 += -2*(m_x[i] - m_a[i])*sum1;
        term3 += sum2*sum2;
        term4 += sum3;
    }

    Ek = -1.0/(2*pow(m_sigma, 4))*(term1 + term2 + term3 + term4) + m_nv/(2*m_sigma*m_sigma);

    return Ek;

}





