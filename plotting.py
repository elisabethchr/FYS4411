import numpy as np
import matplotlib.pyplot as plt
import sys, glob
# from os import listdir
# from os.path import isfile, join

filenames_analytic = glob.glob('data/b/*analytic.txt')
filenames_numeric = glob.glob('data/b/*numeric.txt')

#a1, b2 = np.shape(filenames_analytic)



def convert_to_float(a):
    f = []
    for i in a:
        f.append(float(i))
    return f

def var(filename):
    var = filename.split('_')
    nPart = var[3]
    nDim = var[5]
    nSteps = var[7]
    stepLength = var[9]
    type = var[10].split('.')[0]
    return nPart, nDim, nSteps, stepLength

def variables(filenames):
    nP = []
    nD = []
    nS = []
    sLength = []
    #    t = []
    for file in filenames:
        part, dim, steps, stepLength = var(file)
        nP.append(part)
        nD.append(dim)
        nS.append(steps)
        sLength.append(stepLength)
        #        t.append(type)
    return nP, nD, nS, sLength

nPart_anal, nDim_anal, nSteps_anal, stepLength_anal = variables(filenames_analytic)
nPart_num, nDim_num, nSteps_num, stepLength_num = variables(filenames_numeric)


def filenames(particles, numDimensions, output):
    numParticles_run = len(particles)
    analytic = np.zeros((numParticles_run, numDimensions), dtype=object)
    numeric = np.zeros((numParticles_run, numDimensions), dtype=object)
    part_analytic = np.zeros((numParticles_run, numDimensions), dtype=object)
    part_numeric = np.zeros((numParticles_run, numDimensions), dtype=object)
    dim_analytic = np.zeros((numParticles_run, numDimensions), dtype=object)
    dim_numeric = np.zeros((numParticles_run, numDimensions), dtype=object)

    # organizing filenames according to number of particles and analytic/numeric
    for i in xrange(len(filenames_analytic)):
#        print filenames_analytic[i]
        nPart, nDim, nSteps, steplength = var(filenames_analytic[i])
        if int(nPart)==particles[0]:
            p = 0
        elif int(nPart)==particles[1]:
            p = 1
        elif int(nPart)==particles[2]:
            p = 2
        analytic[p, int(nDim)-1] = filenames_analytic[i]
        part_analytic[p, int(nDim)-1] = nPart

        nPart, nDim, nSteps, steplength = var(filenames_numeric[i])
        if int(nPart)==particles[0]:
            p = 0
        elif int(nPart)==particles[1]:
            p = 1
        elif int(nPart)==particles[2]:
            p = 2
        numeric[p, int(nDim)-1] = filenames_numeric[i]
        part_numeric[p, int(nDim)-1] = nPart
    # print numeric
    # print analytic
    # print part_numeric
    # print part_analytic

    filenames = np.zeros((numParticles_run, numDimensions), dtype=object)
    nParticles = np.zeros((numParticles_run, numDimensions), dtype=object)

    for i in xrange(numParticles_run):
        for j in xrange(numDimensions):
            filenames[i, j] = [numeric[i, j], analytic[i, j]]
            nParticles[i, j] = [part_numeric[i, j], part_analytic[i, j]]

    # plotting
    fig, ax = plt.subplots(numParticles_run, numDimensions)
    for i in xrange(len(particles)):
        for j in xrange(numDimensions):
            data_num = np.loadtxt(filenames[i, j][0], skiprows=1)
            data_anal = np.loadtxt(filenames[i, j][1], skiprows=1)
            alpha_num = data_num[:, 0]
            E_num = data_num[:, 1]
            alpha_anal = data_anal[:, 0]
            E_anal = data_anal[:, 1]
            ax[i, j].plot(alpha_num, E_num, 'r.')
            ax[i, j].plot(alpha_anal, E_anal, 'b--')
            ax[i, j].set_xlabel(r'$\alpha$')
            ax[i, j].set_ylabel(r'$E_{local}$')
            ax[i, j].set_title('nPart = %i, %sD' %(particles[i], j+1))
            ax[i, j].grid(b=True, which='major', color='#666666', linestyle='-')
            # Show the minor grid lines with very faint and almost transparent grey lines
            ax[i, j].minorticks_on()
            ax[i, j].grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
    plt.tight_layout()
#    plt.savefig(output)
    plt.show()

    return filenames, nParticles

particles = [1, 10, 100]
numParticles_run = len(particles)
numDimensions = 3
output = 'data/plots/b/alpha_1_100.pdf'
filenames, nParticles = filenames(particles, 3, output)

print filenames, nParticles
