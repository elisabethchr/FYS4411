#pragma once
#include "../system.h"
#include "optimizer.h"
#include <vector>
#include <armadillo>

using namespace arma;

class StochasticGradientDescent : public Optimizer{
public:
    StochasticGradientDescent(System* system, double eta);

    void computeWeights(arma::vec gradE);

private:
    double m_eta;
    arma::vec m_a;
    arma::vec m_b;
    arma::vec m_w;
};
