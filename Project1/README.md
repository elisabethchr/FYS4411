# Variational Monte Carlo solver

In this project we apply the variational Monte-Carlo method to find the ground state of a dilute hard-sphere bose gas trapped in a harmonic potential. The method is first tested against the closed-form expression for a non-interacting system with a spherically symmetric potential by use of the Metropolis algorithm. It is then expanded to include two-body interactions in an elliptical potential by use of the Metropolis-Hastings algorithm. This is done in one, two and three dimensions for up to 500 particles. Due to correlated measurements, the blocking method is applied in order to present a proper error analysis. The optimal variational parameter for the trial wavefunction is calculated analytically to be $\alpha = 0.5$ for the non-interacting spherically symmetric system. No closed-form solution is obtainable for the elliptical system and so we resort to the steepest descent method to find the value of $\alpha$ that minimizes local energy. Lastly, we plot the one-body density of the interacting and non-interacting system and find that it qualitatively fits our theoretical model.

## Algorithms used:

**Variational Monte Carlo methods**
- Metropolis (system.cpp)
- Metropolis-Hastings (system.cpp)

**Hamiltonians**
- Elliptical harmonic oscillator (Hamiltonians/ellpitical_harmonicoscillator.h/cpp)
- Spherical harmonic oscillator (Hamiltonians/harmonicoscillator.h/cpp)

**Wavefunctions**
- Elliptical gaussian (WaveFunctions/ellipticalguassian.h/cpp)
- Simple gaussian (WaveFunctions/simplegaussian.h/cpp)

**Error analysis**
- Blocking (blocking.py)

**Optimal parameter**
- Gradient descent (main.cpp)
