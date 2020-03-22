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

    double delta = 0.1;
    mat position(m_numberOfParticles, m_numberOfDimensions); position.zeros();

/*
    // Initialize positions
    for(int i=0; i<m_numberOfParticles; i++){
        for(int j=0; j<m_numberOfDimensions; j++){
            position[i, j] = delta*(RandomNumberGenerator(gen) - 0.5);
        }
    }
*/
    double rand;
    // Assign random number with uniform distribution to particles in dimension j
    for (int i=0; i < m_numberOfParticles; i++) {
        std::vector<double> position = std::vector<double>();

        for (int j=0; j < m_numberOfDimensions; j++) {
            /* This is where you should actually place the particles in
             * some positions, according to some rule. Since this class is
             * called random uniform, they should be placed randomly according
             * to a uniform distribution here. However, later you will write
             * more sub-classes of the InitialState class in which the
             * particles are placed in other configurations.
             *
             * Note: For now, the particles are simply placed in positions
             * according to their index in the particles list (this is
             * obviously NOT a good idea).
             */
            rand = delta*(RandomNumberGenerator(gen) - 0.5);


            position.push_back(rand);
        }

        m_particles.push_back(new Particle());
        m_particles.at(i)->setNumberOfDimensions(m_numberOfDimensions);
        m_particles.at(i)->setPosition(position);
    }

}
