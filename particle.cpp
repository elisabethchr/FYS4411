#include "particle.h"
#include <cassert>
#include <armadillo>

using namespace arma;

Particle::Particle() {
}

void Particle::setPosition(const std::vector<double> &position) {
    assert(position.size() == m_numberOfDimensions);        // assert(expression) returns error if position.size() not equal to number of dimensions.
    m_position = position;
}

void Particle::adjustPosition(double change, int dimension) {
    m_position.at(dimension) += change;     //return element at index dimension, e.g. m_position = {1, 5, 2, 4, 3 => m_position.at(3) = 4
}

void Particle::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}
/*
getPosition{position};
    m_getPosition = m_position
    */

