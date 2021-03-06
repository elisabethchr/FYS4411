# FYS4411
This is the course on computational physics: quantum mechanical systems at UiO, taken spring 2020.

## Project 1
# Variational Monte Carlo solver

In this project we apply the variational Monte-Carlo method to find the ground state of a dilute hard-sphere bose gas trapped in a harmonic potential. The method is first tested against the closed-form expression for a non-interacting system with a spherically symmetric potential by use of the Metropolis algorithm. It is then expanded to include two-body interactions in an elliptical potential by use of the Metropolis-Hastings algorithm. This is done in one, two and three dimensions for up to 500 particles. Due to correlated measurements, the blocking method is applied in order to present a proper error analysis. The optimal variational parameter for the trial wavefunction is calculated analytically to be $\alpha = 0.5$ for the non-interacting spherically symmetric system. No closed-form solution is obtainable for the elliptical system and so we resort to the steepest descent method to find the value of $\alpha$ that minimizes local energy. Lastly, we plot the one-body density of the interacting and non-interacting system and find that it qualitatively fits our theoretical model.

## Project 2
# Restricted Boltzmann machine (RBM) applied to the quantum many-body problem

In this project, we employ the restricted Boltzmann machine to find the expected ground-state energy of electrons confined in an isotropic harmonic oscillator potential. The neural-network quan- tum state is optimized by use of stochastic gradient descent. Three different Marcov chain Monte Carlo methods are used to estimate said energy, and their performance is compared. These methods consist of the Metropolis, Metropolis-Hastings and Gibbs sampling procedures.
