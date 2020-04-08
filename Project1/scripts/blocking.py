
# from numpy import log2, zeros, mean, var, sum, loadtxt, arange, \
#                   array, cumsum, dot, transpose, diagonal, floor


import numpy as np
from numpy import log2, zeros, mean, var, sum, loadtxt, arange, array, cumsum, dot, transpose, diagonal, floor
from numpy.linalg import inv
from time import time
import matplotlib.pyplot as plt
import glob, os
# from scipy.interpolate import spline



#filename = '3d_bruteForce_analytic_nParticles_1_nDim_3_nSteps_4194304_dt_0.010000_alpha_0.300000_.txt'  #No dip
#filename = '3d_bruteForce_analytic_nParticles_1_nDim_3_nSteps_16777216_dt_0.010000_alpha_0.300000_.txt'  #No dip
#filename = '3d_bruteForce_analytic_nParticles_1_nDim_3_nSteps_268435456_dt_0.010000_alpha_0.300000_.txt'  #No dip

#filename = 'proper_3d_bruteForce_analytic_nParticles_1_nDim_3_nSteps_1048576_dt_0.010000_alpha_0.300000_steplength0.1.txt'  #No dip
#filename = 'proper_3d_bruteForce_analytic_nParticles_1_nDim_3_nSteps_1048576_dt_0.010000_alpha_0.300000_steplength1.txt'  #No dip

#filename = 'proper_3d_bruteForce_analytic_nParticles_1_nDim_3_nSteps_16777216_dt_0.010000_alpha_0.300000_steplength1.txt'  #No dip
#filename = 'proper_3d_bruteForce_analytic_nParticles_1_nDim_3_nSteps_16777216_dt_0.010000_alpha_0.300000_steplength0.1.txt'  #No dip

#filename = 'proper_3d_bruteForce_analytic_nParticles_1_nDim_3_nSteps_16777216_dt_0.010000_alpha_0.300000_steplength1.txt'
#filename = 'proper_3d_bruteForce_analytic_nParticles_1_nDim_3_nSteps_16777216_dt_0.010000_alpha_0.300000_steplength0.1.txt'
#
#filename = 'proper/1e_AratioNewbruteForce_analytic_nParticles_50_nDim_3_nSteps_131072_dt_0.010000_alpha_0.250000_stepLength_1.000000_.txt'


def type(filename):
    variables = filename.split('_')
    type = variables[3]
    return type

# filename_ideal = 'ideal_onebody_dv2_nPart_1_nDim_3_nSteps_100000000_stepLength_0.100000_nBins_500_non-interacting_analytic.txt'  #No dip
# filename_noJastrow = 'Nonideal_onebody_dv2_nPart_1_nDim_3_nSteps_100000000_stepLength_0.100000_nBins_500_non-interacting_analytic.txt' #No dip
# filename_jastrow = 'plusDR_spherical_onebody_dv2_nPart_20_nDim_3_nSteps_10000000_stepLength_0.100000_nBins_500_interacting_analytic.txt' #No dip


#x = np.loadtxt("3d_bruteForce_analytic_nParticles_1_nDim_3_nSteps_268435456_dt_0.010000_alpha_0.300000_.txt")

x = np.loadtxt('../data/e/1e_importance_analytic_nParticles_1_nDim_3_nSteps_1048576_dt_0.010000_alpha_0.650000_.txt')

def blocking(data1):
    n = len(data1)
    d = int(np.log2(n))
    data = data1
#    mean = np.mean(data)
    sigma_k = np.zeros(d)     # d number of variances (one for each block transformation)
    gamma_k = np.zeros(d)
    e_k0 = np.zeros(d)
    exponent = np.zeros(d)
    EK=0

    # make d number of block transformations
    for i in range(0,d):
        mean = np.mean(data)
        exponent[i] = i
        n = len(data)

        # auto-covariance of energies
       # e_k[i] = auto_covariance(data, mean)

#        if(i==10):
#           EK = auto_covariance(data,mean)
        e_k=0
        # variance of energies
        sigma_k[i] = np.var(data)/n

        # create two arrays containing respecitvely odd and even values of data
        X1 = np.zeros(int(len(data)/2.))
        X2 = np.zeros(int(len(data)/2.))
        for i in range(int(len(data)/2.)):
            X1[i] = data[2*i + 1]
            X2[i] = data[2*i]

        data = 1/2.*(X1 + X2)

#    blockingMean = mean
#    blockingVar = np.average(sigma)
#    #blockingVar = np.sum(sigma) / (float(n-d))
#    blockingErr = np.sqrt(blockingVar)
#    blockingCov = np.average(corr)
    var = sigma_k + e_k

    return var, sigma_k, e_k, exponent, EK

def block(x):
    # preliminaries
    n = len(x); d = int(log2(n)); s, gamma = zeros(d), zeros(d);
    mu = mean(x); t0 = time()

    # estimate the auto-covariance and variances
    # for each blocking transformation
    for i in arange(0,d):
        n = len(x)
        # estimate autocovariance of x
        gamma[i] = (n)**(-1)*sum( (x[0:(n-1)]-mu)*(x[1:n]-mu) )
        # estimate variance of x
        s[i] = np.var(x)
        # perform blocking transformation
        print(np.log2(len(x)))
        x = 0.5*(x[0::2] + x[1::2])

    # generate the test observator M_k ferr2rom the theorem
    M = (cumsum( ((gamma/s)**2*2**arange(1,d+1)[::-1])[::-1] )  )[::-1]

    # we need a list of magic numbers
    q =array([6.634897,9.210340, 11.344867, 13.276704, 15.086272, 16.811894, 18.475307,
              20.090235, 21.665994, 23.209251, 24.724970, 26.216967, 27.688250, 29.141238,
              30.577914, 31.999927, 33.408664, 34.805306, 36.190869, 37.566235, 38.932173,
              40.289360, 41.638398, 42.979820, 44.314105, 45.641683, 46.962942, 48.278236, 49.587884, 50.892181])

    # use magic to determine when we should have stopped blocking
    for k in arange(0,d):
        if(M[k] < q[k]):
            break
    if (k >= d-1):
        print("Warning: Use more data")
    ans = s[k]/2**(d-k)
    print("Runtime: %g sec" % (time()-t0)); print("Blocking Statistics :")
    print("average            iterations      std. error")
    print("%8g %20g %15g" % (mu, k, ans**.5))
    return ans


#var,sigma_k, e_k, exponent, EK = blocking(energy)
filenames = sorted(glob.glob('../data/e/1e_bruteForce_analytic_nParticles_10_nDim_3_nSteps_1048576_dt_0.010000_alpha_*_.txt'), key=os.path.getmtime)

blockingErr = []
for x in filenames:
    print x
    MC = np.loadtxt(x, skiprows=1)
    energy = MC[:,1]

    d = np.log2(len(energy))
    if(d-np.floor(d)!=0):
        print("Number of steps not a power of 2. Truncating measurements.")
        energy = energy[0:2**int(np.floor(d))]
    var = block(energy)
    err = np.sqrt(var)
    blockingErr.append(err)


file = '../data/e/1e_alpha_bruteForce_analytic_nPart_10_nDim_3_stepLength_0.800000_.txt'

header = "Alpha <E> Error Error_b Accepted t_CPU"
#datafile_path = "../data/c/blocking2_ANewRatio_error_50particle_3dim_2^17_steps_bruteForce_stepLength_1.0_"
datafile_path = "../data/e/blocking_error_test_10particle_3dim_2^17_steps_bruteforce_steplength_0.8_"

t = type(file)
#legends[j] = r"$\sigma$ " + t
ofile = datafile_path + t +".txt"

data = np.loadtxt(file, skiprows=1)
alphas = data[:, 0]
energies = data[:, 2]
variance = data[:, 3]
std = data[:, 4]
accepted_ratio = data[:, 5]
t_CPU = data[:, 6]
np.savetxt(ofile, np.column_stack((alphas, energies, std, blockingErr, accepted_ratio, t_CPU)), fmt='%1.2f %1.4f %1.4f %1.4f %1.4f %1.4f', delimiter='   ',  header = header, comments='')

"""
ans = block(energy)
a = np.empty(len(var))
a.fill(ans)

mean = np.mean(energy)
varC =np.zeros(len(sigma_k))
#EK = auto_covariance(energy_a, mean)
for l in range(0,int(d)-1):
    varC[l] = sigma_k[10] + EK


print(var)
print(exponent)


fig0 = plt.figure(0)
ax = fig0.gca()
# legends.append(r'$\sigma_b$ numeric')
# legends.append(r'$\sigma_b$ analytic')


#plt.plot(exponent, sigma_k, linestyle = '-',label='$\sigma_k^2/n_k$')
plt.plot(exponent, var,'r',linestyle = '--',label='$\sigma_k^2/n_k$')
plt.plot(exponent,a,'b',linestyle ='--', label = '$\sigma_b^2$')
#plt.plot(exponent,varC,linestyle='--',label='$Var(\\overline{X}$')
plt.xlabel('Blocking iteration, $k$', size=12)
plt.ylabel('Variance', size=12)



plt.legend()
ax.set_xticks(np.arange(0, int(np.floor(d)), 1))

plt.grid(b=True, which='major', color='#666666', linestyle='-')
# Show the minor grid lines with very faint and almost transparent grey lines
#plt.minorticks_on()
#plt.grid(b=true, which='minor', color='#999999', linestyle='-', alpha=0)
#plt.tight_layout()
#plt.grid('on')
# plt.savefig('data/plots/d/blocking_error_1particle_3dim_2^19_steps_bruteForce.pdf')
plt.show()
"""
