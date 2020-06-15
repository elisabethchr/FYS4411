# Restricted Boltzmann machine (RBM) applied to the quantum many-body problem

In this project, we employ the restricted Boltzmann machine to find the expected ground-state energy of electrons confined in an isotropic harmonic oscillator potential. The neural-network quan- tum state is optimized by use of stochastic gradient descent. Three different Marcov chain Monte Carlo methods are used to estimate said energy, and their performance is compared. These methods consist of the Metropolis, Metropolis-Hastings and Gibbs sampling procedures.


## Algorithms used:

**Monte Carlo solvers**
- Metropolis (system.cpp)
- Metropolis-Hastings (system.cpp)
- Gibbs sampling (system.cpp)

**Hamiltonians**
- Spherical harmonic oscillator (Hamiltonians/harmonicoscillator.h/cpp)

**Wavefunctions**
- Neural quantum state (WaveFunctions/neuralquantumstate.h/cpp)

**Error analysis**
- Blocking (blocking.py)

**Optimizing weights**
- Stochastic gradient descent (Optimizers/stochasticgradientdescent.h/cpp)

