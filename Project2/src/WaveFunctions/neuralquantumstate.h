#pragma once
#include <armadillo>
#include "wavefunction.h"

class NeuralQuantumState : public WaveFunction {
public:
    NeuralQuantumState(class System* system, int n_hidden, int n_visible, int part, int dim, double sigma);
    double evaluate(arma::vec position);

    arma::vec set_X(arma::vec X){ m_x = X; }
    arma::vec set_a(arma::vec a){ m_a = a; }
    arma::vec set_b(arma::vec b){ m_b = b; }
    arma::mat set_w(arma::mat w) {m_w = w; }

    arma::vec get_X(){ return m_x; }
    arma::vec get_a(){ return m_a; }
    arma::vec get_b(){ return m_b; }
    arma::mat get_w(){ return m_w; }

private:
    void setupInitialState();
    void setupWeights(){}
    void setupPositions(){}

    arma::vec O;        // vector exponent in wavefunction
    arma::vec m_x;      // visible nodes (i.e. position)
    arma::vec m_a;      // visible bias
    arma::vec m_b;      // hidden bias
    arma::mat m_w;      // interaction of biases

    double m_term1;
    double m_term2;
    double m_sigma;

};


