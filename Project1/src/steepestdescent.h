#pragma once
#include <ctime>
#include <armadillo>

class SteepestDescent{
public:
    SteepestDescent(class System* system);

private:
    class System* m_system = nullptr;
}
