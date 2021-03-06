#pragma once
#include <armadillo>
#include "wavefunction.h"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(class System* system, std::vector<double> alpha);
    double evaluate(std::vector<class Particle*> particles);
    arma::vec computeGradient(std::vector<class Particle*> particles, int i);
    double computeDoubleDerivative_numeric(std::vector<class Particle*> particles);
    double computeDoubleDerivative_analytic(std::vector<class Particle *> particles);
    double getDistance(std::vector<class Particle *> particles, int particle_i, int particle_j){}
    double computeLaplacian(std::vector<class Particle*> particles){}
    double computeDerivativePsi_alpha(std::vector<class Particle*> particles){}

protected:
    double SingleParticleFunction(std::vector<class Particle *> particles, int particle){}
    double correlationWaveFunction(std::vector<class Particle *> particles, double distance){}
    double computeLaplacian_SPF(int particle_i){}
    arma::vec computeGradient_SPF(int particle_i){}
    double computeLaplacian_correlation(std::vector<class Particle *> particles, int particle_i){}
    arma::vec computeGradient_correlation(std::vector<class Particle *> particles, int particle_i){}
    double computeDerivative_correlation(int particle_i, int particle_j){}
    double computeDoubleDerivative_correlation(int particle_i, int particle_j){}

};
