TEMPLATE = app
CONFIG  += console c++11
CONFIG  -= app_bundle
CONFIG  -= qt

LIBS += -llapack -lblas -larmadillo

SOURCES += main.cpp \
    WaveFunctions/neuralquantumstate.cpp \
    system.cpp \
    Hamiltonians/hamiltonian.cpp \
    Hamiltonians/harmonicoscillator.cpp \
    particle.cpp \
    WaveFunctions/wavefunction.cpp \
    InitialStates/initialstate.cpp \
    InitialStates/initalizestate.cpp \
    Math/random.cpp \
    sampler.cpp \

HEADERS += \
    WaveFunctions/neuralquantumstate.h \
    system.h \
    Hamiltonians/hamiltonian.h \
    Hamiltonians/harmonicoscillator.h \
    particle.h \
    WaveFunctions/wavefunction.h \
    InitialStates/initialstate.h \
    InitialStates/initalizestate.h \
    Math/random.h \
    sampler.h \
