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


def filenames(particles, numDimensions):
    numParticles_run = len(particles)
    analytic = np.zeros((numParticles_run, numDimensions), dtype=object)
    numeric = np.zeros((numParticles_run, numDimensions), dtype=object)
    part_analytic = np.zeros((numParticles_run, numDimensions), dtype=object)
    part_numeric = np.zeros((numParticles_run, numDimensions), dtype=object)
    dim_analytic = np.zeros((numParticles_run, numDimensions), dtype=object)
    dim_numeric = np.zeros((numParticles_run, numDimensions), dtype=object)

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
    print numeric
    print analytic
    print part_numeric
    print part_analytic

    filenames = np.zeros((numParticles_run, numDimensions), dtype=object)
    nParticles = np.zeros((numParticles_run, numDimensions), dtype=object)
#    nDim = np.zeros()

    for i in xrange(numParticles_run):
        for j in xrange(numDimensions):
            filenames[i, j] = [numeric[i, j], analytic[i, j]]
            nParticles[i, j] = [part_numeric[i, j], part_analytic[i, j]]
    return filenames, nParticles

particles = [1, 10]
numParticles_run = len(particles)
numDimensions = 3

filenames, nParticles = filenames(particles, 3)

print filenames, nParticles

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
#        ax[i, j].set_title('nPart = %s' %nParticles[i, j][0])
#        ax[i, j].grid('on')
        ax[i, j].grid(b=True, which='major', color='#666666', linestyle='-')
        ax[i, j].minorticks_on()
        ax[i, j].grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)

# Show the minor grid lines with very faint and almost transparent grey lines
plt.tight_layout()
plt.savefig('data/plots/b/alpha_1_100.pdf')
plt.show()

# """
# filename1 = 'data/1b_alpha_nPart_1_nDim_1_nSteps_100000_stepLength_1.000000_numeric.txt'
# filename2 = "data/1b_alpha_nPart_1_nDim_1_nSteps_100000_stepLength_1.000000_analytic.txt"
# filename3 = 'data/1b_alpha_nPart_1_nDim_1_nSteps_100000_stepLength_0.500000_numeric.txt'
# filename4 = "data/1b_alpha_nPart_1_nDim_1_nSteps_100000_stepLength_0.500000_analytic.txt"
# filename5 = 'data/1b_alpha_nPart_1_nDim_1_nSteps_100000_stepLength_0.100000_numeric.txt'
# filename6 = "data/1b_alpha_nPart_1_nDim_1_nSteps_100000_stepLength_0.100000_analytic.txt"
#
# nPart, nDim, nSteps, stepLength1 = var(filename1)
# nPart, nDim, nSteps, stepLength2 = var(filename3)
# nPart, nDim, nSteps, stepLength3 = var(filename5)
#
# data_num1 = np.loadtxt(filename1, skiprows=1)
# data_anal1 = np.loadtxt(filename2, skiprows=1)
# data_num05 = np.loadtxt(filename3, skiprows=1)
# data_anal05 = np.loadtxt(filename4, skiprows=1)
# data_num01 = np.loadtxt(filename4, skiprows=1)
# data_anal01 = np.loadtxt(filename5, skiprows=1)
#
# alpha_num1 = data_num1[:, 0]
# E_num1 = data_num1[:, 1]
# alpha_anal1 = data_anal1[:, 0]
# E_anal1 = data_anal1[:, 1]
# alpha_num05 = data_num05[:, 0]
# E_num05 = data_num05[:, 1]
# alpha_anal05 = data_anal05[:, 0]
# E_anal05 = data_anal05[:, 1]
# alpha_num01 = data_num01[:, 0]
# E_num01 = data_num01[:, 1]
# alpha_anal01 = data_anal01[:, 0]
# E_anal01 = data_anal01[:, 1]
#
# plt.plot(alpha_num1, E_num1, 'r.')
# plt.plot(alpha_anal1, E_anal1, 'b--')
# plt.plot(alpha_num05, E_num05, 'g+')
# plt.plot(alpha_anal05, E_anal05, 'y--')
# plt.plot(alpha_num01, E_num01, 'k^')
# plt.plot(alpha_anal01, E_anal01, 'm--')
# plt.legend(["Numeric %s" %stepLength1, "Analytic %s" %stepLength1, "Numeric %s" %stepLength2, "Analytic %s" %stepLength2, "Numeric %s" %stepLength3 , "Analytic %s" %stepLength3])
# plt.xlabel(r'$\alpha$')
# plt.ylabel(r'$E_L$')
# plt.title("Steps = %s" %nSteps)
# plt.tight_layout()
# plt.grid('on')
# plt.savefig("data/plots/stepLengths.pdf")
# plt.show()
# """
