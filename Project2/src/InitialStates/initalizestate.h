#pragma once
#include "initialstate.h"

class InitializeState : public InitialState {
public:
    InitializeState(System* system, int n_hidden, int n_visible, int dim);
    void setupInitialState();
    void setupPositions();
    void setupWeights(){}
};

