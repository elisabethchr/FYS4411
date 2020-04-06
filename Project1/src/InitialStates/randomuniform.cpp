#include "randomuniform.h"
#include <iostream>
#include <cassert>
#include <armadillo>
#include "Math/random.h"
#include "../particle.h"
#include "../system.h"

using std::cout;
using std::endl;
using namespace std;
using namespace arma;

RandomUniform::RandomUniform(System* system, int numberOfDimensions, int numberOfParticles)  : InitialState(system) {
    assert(numberOfDimensions > 0 && numberOfParticles > 0);
    m_numberOfDimensions = numberOfDimensions;
    m_numberOfParticles  = numberOfParticles;

    /* The Initial State class is in charge of everything to do with the
     * initialization of the system; this includes determining the number of
     * particles and the number of dimensions used. To make sure everything
     * works as intended, this information is passed to the system here.
     */
    m_system->setNumberOfDimensions(numberOfDimensions);
    m_system->setNumberOfParticles(numberOfParticles);
    setupInitialState();
}

void RandomUniform::setupInitialState() {

    // Random number generator
    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);

    double delta = 0.5;
    mat position(m_numberOfParticles, m_numberOfDimensions); position.zeros();

    double rand;
    // Assign random number with uniform distribution to particles in dimension j
    for (int i=0; i < m_numberOfParticles; i++) {
        std::vector<double> position = std::vector<double>();

        for (int j=0; j < m_numberOfDimensions; j++) {
            rand = delta*(RandomNumberGenerator(gen) - 0.5);
            position.push_back(rand);
        }

        m_particles.push_back(new Particle());
        m_particles.at(i)->setNumberOfDimensions(m_numberOfDimensions);
        m_particles.at(i)->setPosition(position);
    }

}
