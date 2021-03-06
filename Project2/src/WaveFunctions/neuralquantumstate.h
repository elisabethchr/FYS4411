#pragma once
#include <armadillo>
#include "../system.h"
#include "wavefunction.h"

class NeuralQuantumState : public WaveFunction {
public:
    NeuralQuantumState(class System* system, int n_hidden, int n_visible, int part, int dim, double sigma, bool gaussian, double initialization);
    double evaluate(arma::vec position);

    arma::vec set_X(arma::vec X){ m_x = X; }
    arma::vec set_h(arma::vec h){ m_h = h; }
    arma::vec set_a(arma::vec a){ m_a = a; }
    arma::vec set_b(arma::vec b){ m_b = b; }
    arma::mat set_w(arma::mat w) {m_w = w; }

    arma::vec get_X(){ return m_x; }
    arma::vec get_h(){ return m_h; }
    arma::vec get_a(){ return m_a; }
    arma::vec get_b(){ return m_b; }
    arma::mat get_w(){ return m_w; }

    double computeDerivative_analytic(arma::vec position, int coor);
    double computeDoubleDerivative_analytic();
    double getDistance(int p, int q);
    double getSigma(){ return m_sigma; }
    double getInitializationInterval() {return m_initialization; }
    double getNumberVisibleNodes(){ return m_nv; }
    double getNumberHiddenNodes(){ return m_nh; }
    bool getGaussian(){ return m_gaussian; }

private:
    void setupInitialState();

    double sigmoid(double x); // the logistic function
    double v(int j);          // for simplification
    arma::mat hadamardProd(arma::mat w); // calculates the hadamard product of w with itself

    arma::vec O;        // vector exponent in wavefunction
    arma::vec m_x;      // visible nodes (i.e. position)
    arma::vec m_h;      // hidden nodes
    arma::vec m_a;      // visible bias
    arma::vec m_b;      // hidden bias
    arma::mat m_w;      // interaction of biases

    double m_term1;
    double m_term2;
    double m_sigma;
    double m_nv;
    double m_nh;
    double m_initialization;

    bool m_gaussian;      // determines distribution of weights (either uniform or normal)
};
