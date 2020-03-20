# Calculate the blocking method on different solvers, write calculations to file and plot errors

import numpy as np
import matplotlib.pyplot as plt
import sys, glob, os

def alpha(filename):
    variables = filename.split('_')
    alpha = variables[12]
    return float(alpha)

def type(filename):
    variables = filename.split('_')
    type = variables[3]
    return type

def dimension(filename):
    variables = filename.split('_')
    dim = variables[7]
    return int(dim)

def auto_covariance(data, mu):
    n = len(data)
    autoCov = np.sum((data[0:n-1] - mu)*(data[1:n] - mu))
    return autoCov/float(n)

# blocking method
def blocking(data1):
    n = len(data1[:, 0])
    d = int(np.log2(n))
    data = data1[:, 1]
    mean = np.mean(data)
    sigma = np.zeros(d)     # d number of variances (one for each block transformation)
    gamma_k = np.zeros(d)
    corr = np.zeros(d)

    # make d number of block transformations
    for i in xrange(d-1):
        n = len(data)

        # auto-covariance of energies
        corr[i] = auto_covariance(data, mean)

        # variance of energies
        sigma[i] = np.var(data)

        # create two arrays containing respecitvely odd and even values of data
        X1 = np.zeros(int(len(data)/2.))
        X2 = np.zeros(int(len(data)/2.))
        for i in xrange(int(len(data)/2.)):
            X1[i] = data[2*i + 1]
            X2[i] = data[2*i]

        data = 1/2.*(X1 + X2)

    blockingMean = mean
    blockingVar = np.average(sigma)
    #blockingVar = np.sum(sigma) / (float(n-d))
    blockingErr = np.sqrt(blockingVar)
    blockingCov = np.average(corr)

    return blockingMean, blockingVar, blockingCov, blockingErr


# get filenames and sorted according to time created
timesteps_numeric = sorted(glob.glob('data/d/1d_bruteForce_numeric_nParticles_1_nDim_3*'), key=os.path.getmtime)
timesteps_analytic = sorted(glob.glob('data/d/1d_bruteForce_analytic_nParticles_1_nDim_3*'), key=os.path.getmtime)

filenames = [timesteps_numeric, timesteps_analytic]

blocking_alphas = [[], []]
blockingMean = [[], []]
blockingVar = [[], []]
blockingCov = [[], []]
blockingErr = [[], []]

# variance from blocking, i.e. correlated standard deviations
for i in xrange(len(filenames)):
    for j in xrange(len(filenames[i])):
        file = filenames[i][j]
        print file
        data = np.loadtxt(file, skiprows=1)

        a = alpha(file)
        blocking_alphas[i].append(a)

        mean, var, cov, err = blocking(data)
        blockingMean[i].append(mean)
        blockingVar[i].append(var)
        blockingErr[i].append(err)
        blockingCov[i].append(cov)

print blockingErr

# variance from data, i.e. uncorrelated standard deviation
dim = 3
filenames = glob.glob('data/d/1d_alpha_bruteForce*')

filenames_alpha = []
for i in filenames:
    if dimension(i) == dim:
        filenames_alpha.append(i)

print filenames_alpha
legends = [[]]*len(filenames_alpha)
alphas = [[]]*len(filenames_alpha)
energies = [[]]*len(filenames_alpha)
variance = [[]]*len(filenames_alpha)
std = [[]]*len(filenames_alpha)
legends = [[]]*len(filenames_alpha)
accepted_ratio = [[]]*len(filenames_alpha)
t_CPU = [[]]*len(filenames_alpha)

header = "Alpha <E> Error Error_b Accepted t_CPU"
datafile_path = "data/d/blocking_error_1particle_3dim_2^19_steps_bruteForce_"

j=0
for i in filenames_alpha:
    print i
    t = type(i)
    legends[j] = r"$\sigma$ " + t
    ofile = datafile_path + t +".txt"

    data = np.loadtxt(i, skiprows=1)

    alphas[j] = data[:, 0]
    energies[j] = data[:, 2]
    variance[j] = data[:, 3]
    std[j] = data[:, 4]
    accepted_ratio[j] = data[:, 5]
    t_CPU[j] = data[:, 6]

    np.savetxt(ofile, np.column_stack((alphas[j], energies[j], std[j], blockingErr[j], accepted_ratio[j], t_CPU[j])), fmt='%1.2f %1.4f %1.4f %1.4f %1.4f %1.4f', delimiter='   ',  header = header, comments='')

    plt.plot(alphas[j], std[j], linestyle = '--', marker = '+')
    j+=1

legends.append(r'$\sigma_b$ numeric')
legends.append(r'$\sigma_b$ analytic')
plt.plot(blocking_alphas[0], blockingErr[0], linestyle = '--', marker = '+')
plt.plot(blocking_alphas[1], blockingErr[1], linestyle = '--', marker = '+')
plt.title('Brute force')
plt.xlabel(r'$\alpha$', size=12)
plt.ylabel(r'$\sigma$', size=12)
plt.legend(legends)
plt.grid(b=True, which='major', color='#666666', linestyle='-')
# Show the minor grid lines with very faint and almost transparent grey lines
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.tight_layout()
plt.grid('on')
plt.savefig('data/plots/d/blocking_error_1particle_3dim_2^19_steps_bruteForce.pdf')
plt.show()
