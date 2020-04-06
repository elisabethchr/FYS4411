#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <iomanip>
#include <cstdlib>
#include <algorithm>
#include "conjugategradient.h"
#include "system.h"

int nIterations = 50;
double Energy, EnergyDerivative, alpha_guess, alphagradient, eta, alpha_prev;

eta = 0.;
alpha_guess = 0.9;
alpha.push_back(alpha_guess);

System* system = new System();
system->setCalculation              (numeric);
system->setTimeSteps                (timestep);
system->setSolver                   (bruteForce);

    for (int iter=0; iter<nIterations+1; iter++){
        alpha_prev = alpha[0];
        cout << "*************************************************************" << endl;
        cout << "|" << endl;
        cout << "|" << endl;
        cout << "|" << endl;
        cout << "| Iteration no. " << iter << ", alpha = " << alpha[0] << endl;
        cout << "|" << endl;
        cout << "|" << endl;
        cout << "|" << endl;
        cout << "*************************************************************" << endl;

        system->setHamiltonian              (new EllipticalHarmonicOscillator(system, omega));
        system->setWaveFunction             (new EllipticalGaussian(system, alpha, beta, a));
    //    system->setHamiltonian              (new HarmonicOscillator(system, omega));
    //    system->setWaveFunction             (new SimpleGaussian(system, alpha));
        system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
        system->setEquilibrationFraction    (equilibration);
        system->setStepLength               (stepLength);
        system->runMetropolisSteps          (MC_cycles);

        Energy = system->getEnergy();
        EnergyDerivative = system->getEnergyDerivative();
        alphagradient = EnergyDerivative;
        alpha[0] = alpha_prev - eta*alphagradient;


    }
