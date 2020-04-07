#pragma once
#include <vector>

class ConjugateGradient {
public:
    ConjugateGradient(class System* system);
    void SteepestDescent();

protected:
    class System* m_system = nullptr;

};

