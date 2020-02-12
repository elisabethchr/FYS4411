#pragma once
#include "wavefunction.h"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(class System* system, double alpha);
    double evaluate(std::vector<class Particle*> particles, double alpha);
    double computeDoubleDerivative(std::vector<class Particle*> particles);
};
