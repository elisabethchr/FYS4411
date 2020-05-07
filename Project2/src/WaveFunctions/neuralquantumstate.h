#pragma once
#include <armadillo>
#include "wavefunction.h"

class NeuralQuantumState : public WaveFunction {
public:
    NeuralQuantumState(class System* system);
    double evaluate(std::vector<class Particle*> particles);

    arma::vec computeGradient(std::vector<class Particle*> particles, int i){}
    double computeLaplacian(std::vector<class Particle*> particles){}
    double computeDoubleDerivative_numeric(std::vector<class Particle*> particles){}
    double getDistance(std::vector<class Particle *> particles, int particle_i, int particle_j){}
    double computeDerivativePsi_alpha(std::vector<class Particle*> particles){}

};
