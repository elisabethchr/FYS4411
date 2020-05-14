#pragma once
#include "initialstate.h"

class InitializeState : public InitialState {
public:
    InitializeState(System* system, int numberOfDimensions, int numberOfParticles);
    void setupInitialState();
};

