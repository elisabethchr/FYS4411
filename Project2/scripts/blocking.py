import numpy as np
from numpy import log2, zeros, mean, var, sum, loadtxt, arange, array, cumsum, dot, transpose, diagonal, floor
from numpy.linalg import inv
from time import time
import matplotlib.pyplot as plt
import glob, os
# from scipy.interpolate import spline

def block(x):
    # preliminaries
    n = len(x); d = int(log2(n)); s, gamma = zeros(d), zeros(d);
    mu = mean(x); t0 = time()
    print d

    # estimate the auto-covariance and variances
    # for each blocking transformation
    for i in arange(0,d):
        print "Transformation", i
        n = len(x)
        print "length of block:", n
        # estimate autocovariance of x
        gamma[i] = (n)**(-1)*sum( (x[0:(n-1)]-mu)*(x[1:n]-mu) )
        # estimate variance of x
        s[i] = np.var(x)
        # perform blocking transformation
        print np.log2(len(x))
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
