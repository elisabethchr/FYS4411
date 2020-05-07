#pragma once
#include <vector>
#include <random>

using namespace std;

class InitialState {
public:
    InitialState(class System* system);
    std::vector<class Particle*> getParticles() { return m_particles; }

    virtual void setupInitialState() = 0;
    virtual void setupPositions() = 0;
    virtual void setupWeights() = 0;

protected:
    class System* m_system = nullptr;
    std::vector<Particle*> m_particles;// = std::vector<Particle*>();
    int m_nh = 0;
    int m_nv = 0;
    int m_dim = 0;
    std::mt19937_64 m_randomEngine; // For the distributions
};

