#pragma once
#include <vector>
#include <armadillo>
//#include "../system.h"

using namespace arma;

class Optimizer {
public:
    Optimizer(class System* system);

    virtual void computeWeights(arma::vec gradE) = 0;

protected:
    class System* m_system = nullptr;
};
