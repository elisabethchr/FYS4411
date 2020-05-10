#include "../system.h"
#include "optimizer.h"
#include <vector>

class StochasticGradientDescent : public Optimizer{
public:
    StochasticGradientDescent(System* system, double eta);

    void computeWeights(arma::vec gradE, int cycle);

private:
    double m_eta;
    arma::vec m_a;
    arma::vec m_b;
    arma::vec m_w;
};
