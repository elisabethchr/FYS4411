TEMPLATE = app
CONFIG  += console c++11
CONFIG  -= app_bundle
CONFIG  -= qt

SOURCES += main.cpp \
    system.cpp \
    Hamiltonians/hamiltonian.cpp \
    Hamiltonians/harmonicoscillator.cpp \
    Hamiltonians/elliptical_harmonicoscillator.cpp \
    particle.cpp \
    WaveFunctions/wavefunction.cpp \
    InitialStates/initialstate.cpp \
    InitialStates/randomuniform.cpp \
    Math/random.cpp \
    sampler.cpp \
    WaveFunctions/simplegaussian.cpp \
    WaveFunctions/ellipticalgaussian.cpp \

HEADERS += \
    system.h \
    Hamiltonians/hamiltonian.h \
    Hamiltonians/harmonicoscillator.h \
    Hamiltonians/elliptical_harmonicoscillator.h \
    particle.h \
    WaveFunctions/wavefunction.h \
    InitialStates/initialstate.h \
    InitialStates/randomuniform.h \
    Math/random.h \
    sampler.h \
    WaveFunctions/simplegaussian.h \
    WaveFunctions/ellipticalgaussian.h \




#TEMPLATE = app
#CONFIG += console c++11
#CONFIG -= app_bundle
#CONFIG -= qt

#SOURCES +=
##    main.cpp \
##    Project1/Hamiltonians/hamiltonian.cpp \
##    Project1/Hamiltonians/harmonicoscillator.cpp \
##    Project1/InitialStates/initialstate.cpp \
##    Project1/InitialStates/randomuniform.cpp \
##    Project1/Math/random.cpp \
##    Project1/WaveFunctions/simplegaussian.cpp \
##    Project1/WaveFunctions/wavefunction.cpp \
##    Project1/main.cpp \
##    Project1/particle.cpp \
##    Project1/sampler.cpp \
##    ../../FYS4411/Project1/system.cpp \


#    Hamiltonians/harmonicoscillator.cpp \
#    InitialStates/initialstate.cpp \
#    InitialStates/randomuniform.cpp \
#    Math/random.cpp \
#    WaveFunctions/simplegaussian.cpp \
#    WaveFunctions/wavefunction.cpp \
#    main.cpp \
#    particle.cpp \
#    sampler.cpp \
#    system.cpp

#HEADERS += \

##    Project1/Hamiltonians/hamiltonian.h \
##    Project1/Hamiltonians/harmonicoscillator.h \
##    Project1/InitialStates/initialstate.h \
##    Project1/InitialStates/randomuniform.h \
##    Project1/Math/random.h \
##    Project1/WaveFunctions/simplegaussian.h \
##    Project1/WaveFunctions/wavefunction.h \
##    Project1/particle.h \
##    Project1/sampler.h \
##    Project1/system.h \

#    Hamiltonians/hamiltonian.h \
#    Hamiltonians/harmonicoscillator.h \
#    InitialStates/initialstate.h \
#    InitialStates/randomuniform.h \
#    Math/random.h \
#    WaveFunctions/simplegaussian.h \
#    WaveFunctions/wavefunction.h \
#    particle.h \
#    sampler.h \
#    system.h
