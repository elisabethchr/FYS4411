#pragma once
#include <armadillo>
#include "wavefunction.h"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(class System* system, std::vector<double> alpha);
//    SimpleGaussian(class System* system, double alpha);
    double evaluate(std::vector<class Particle*> particles);
    arma::vec computeGradient(std::vector<class Particle*> particles, int i);
    double computeDoubleDerivative_numeric(std::vector<class Particle*> particles);
    double computeDoubleDerivative_analytic(std::vector<class Particle *> particles);
};
