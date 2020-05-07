#include <cmath>
#include <cassert>
#include <armadillo>
#include <random>
#include "neuralquantumstate.h"
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

using namespace std;
using namespace arma;

NeuralQuantumState::NeuralQuantumState(System* system) :
    WaveFunction(system) {
}

double NeuralQuantumState::evaluate(std::vector<class Particle*> particles){

}
