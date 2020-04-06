import numpy as np
import glob, sys, os

timesteps_analytic = sorted(glob.glob('../data/e/1e_bruteForce_analytic_nParticles_10*'), key=os.path.getmtime)
filenames = sorted(glob.glob('../data/e/1e_alpha_bruteForce_analytic_nPart_*_nDim_3_stepLength_0.100000_*'))

print filenames #okay
