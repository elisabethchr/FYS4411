#pragma once
#include <vector>

class Hamiltonian {
public:
    Hamiltonian(class System* system);
    virtual double computeLocalEnergy(std::vector<class Particle*> particles) = 0;
    virtual double computeLocalEnergy_analytic(std::vector<class Particle*> particles) = 0;
    //Remove later: This is a pure virtual function, meaning any derived classes e.g. HarmonicOscillator MUST overwrite it.
    //In a normal virtual function one would simply omit the "=0", and any derived classes would use the
    //implementation of the base class if it has not been overwritten. Pure virtual functions are good for
    //making sure every derived class has its own implementation of the function.

protected:
    class System* m_system = nullptr;
};

