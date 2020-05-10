TEMPLATE = app
CONFIG  += console c++11
CONFIG  -= app_bundle
CONFIG  -= qt

SOURCES += main.cpp \
    Optimizers/optimizer.cpp \
    Optimizers/stochasticgradientdescent.cpp \
    system.cpp \
    particle.cpp \
    sampler.cpp \
    Hamiltonians/hamiltonian.cpp \
    Hamiltonians/harmonicoscillator.cpp \
    WaveFunctions/wavefunction.cpp \
    WaveFunctions/neuralquantumstate.cpp \
    InitialStates/initialstate.cpp \
    InitialStates/initalizestate.cpp \
    Math/random.cpp \

HEADERS += \
    Optimizers/optimizer.h \
    Optimizers/stochasticgradientdescent.h \
    system.h \
    particle.h \
    sampler.h \
    Hamiltonians/hamiltonian.h \
    Hamiltonians/harmonicoscillator.h \
    WaveFunctions/wavefunction.h \
    WaveFunctions/neuralquantumstate.h \
    InitialStates/initialstate.h \
    InitialStates/initalizestate.h \
    Math/random.h \
