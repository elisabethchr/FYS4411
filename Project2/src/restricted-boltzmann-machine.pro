TEMPLATE = app
CONFIG  += console c++11
CONFIG  -= app_bundle
CONFIG  -= qt

LIBS += -llapack -lblas -larmadillo

SOURCES += main.cpp \
    Optimizers/optimizer.cpp \
    Optimizers/stochasticgradientdescent.cpp \
    system.cpp \
    sampler.cpp \
    Hamiltonians/hamiltonian.cpp \
    Hamiltonians/harmonicoscillator.cpp \
    WaveFunctions/wavefunction.cpp \
    WaveFunctions/neuralquantumstate.cpp \
    Math/random.cpp \

HEADERS += \
    Optimizers/optimizer.h \
    Optimizers/stochasticgradientdescent.h \
    system.h \
    sampler.h \
    Hamiltonians/hamiltonian.h \
    Hamiltonians/harmonicoscillator.h \
    WaveFunctions/wavefunction.h \
    WaveFunctions/neuralquantumstate.h \
    Math/random.h \
