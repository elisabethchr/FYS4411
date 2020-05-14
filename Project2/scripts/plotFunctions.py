import numpy as np
import matplotlib.pyplot as plt
from blocking import block
import glob, os
from decimal import Decimal
import csv

#********************************************************
# Functions for obtaining given variables for a dataset
#********************************************************

# example filename:
# filename = 2b_RBMcycles_bruteForce_x_nv_x_nh_x_nMCsteps_x_nPart_x_nDim_x_stepL_x_eta_x_RBMcycles_x.txt
# where x represents an integer or float number

def distribution(filename):
    dist = filename.split('_')[3]
    return dist

def solver(filename):
    solv = filename.split('_')[2]
    return solv

def nNodes(filename):
    a = filename.split('_')
    nVisible = a[5]
    nHidden = a[7]
    return int(nVisible), int(nHidden)

def nMCsteps(filename):
    nSteps = filename.split('_')[9]
    return int(nSteps)

def nPart(filename):
    part = filename.split('_')[11]
    return int(part)

def nDim(filename):
    dim = filename.split('_')[13]
    return int(dim)

def stepLength(filename):
    stepL = filename.split('_')[15]
    return float(stepL)

def learningRate(filename):
    rate = filename.split('_')[17]
    return float(rate)

def RBMCycles(filename):
    cycles = filename.split('_')[19]
    return int(cycles)

"""
**********************************************************

Plotting of local energies as a function of RBM cycles

**********************************************************
"""

#***************************************************
# Uniform vs normal weight distribution
#***************************************************

# Function for plotting energies vs. RBM cycles for distinct learning rate, with distribution as legend
def RBM_distribution(filenames, datafile_path, ofile):
    distributions = []
    mean = []
    variance = []
    error = []
    rates = []
    colors = ['b', 'r', 'k', 'm', 'g']
    c=0
    ofile = ofile+str(learningRate(filenames[0]))
    for file in filenames:
#        print "File:",file
        data = np.loadtxt(file, skiprows=1)
        cycles = data[1:, 0]
        energies = data[1:, 1]
        dist = distribution(file)
        mu = np.mean(energies)
        sigma2 = np.var(energies, ddof = 1)   # sample variance (assuming the energies from file is only a sample of all energies we could have been looking at, i.e. more RBM cycles)
        eta = learningRate(file)
        distributions.append(dist)
        mean.append(mu)
        variance.append(sigma2)
        error.append(np.sqrt(sigma2))
        plt.plot(cycles, energies, linestyle='--', marker = '.', color=colors[c], label=dist)
        c+=1
    plt.xlabel('cycles')
    plt.ylabel(r'$E_L$')
    plt.title(r'$E_L$ dependence on distributions of weights, $\eta = %.1E$' %Decimal(eta))
    plt.legend()
    # Show the minor grid lines with very faint and almost transparent grey lines
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
    plt.tight_layout()
    plt.savefig(datafile_path+ofile+'_.pdf')
    plt.show()

    header = "Distribution  <E_L>   sigma2  error"

    # write to file
    dataset = np.array([distributions, mean, variance, error])
    dataset = dataset.T
    print dataset

    datafile = datafile_path+ofile+'_.txt'
    with open(datafile, "w+") as datafile:
        np.savetxt(datafile, dataset, fmt = ['%s', '%s', '%s', '%s'], header=header, delimiter = '\t')

    return dataset


#****************************************************
# Energies vs. RBM cycles, variables as legends
# (learning rate and number of hidden nodes)
#****************************************************

# Function for plotting energies vs. RBM cycles for distinct distributions, with learning rate as legend
def plotRBM_learningrate(filenames, datafile_path, ofile):
    dist = distribution(filenames[0])
    if (dist=="uniform"):
        linestyle = '--'
    elif (dist=="normal"):
        linestyle = '-'
    for file in filenames:
        print "File:", file
        data = np.loadtxt(file, skiprows=1)
        cycles = data[1:, 0]
        energies = data[1:, 1]
        eta = learningRate(file)
        plt.plot(cycles, energies, linestyle=linestyle, label="%.1E" %Decimal(eta))
    plt.xlabel("cycles")
    plt.ylabel(r'$E_L$')
    plt.title("Energy dependece on learning rates")
    plt.legend()
    # Show the minor grid lines with very faint and almost transparent grey lines
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
    plt.tight_layout()
    plt.savefig(datafile_path+ofile+"_"+dist+'_.pdf')
    plt.show()

# Function for plotting energies vs. RBM cycles for distinct distributions, with hidden nodes as legend
def plotRBM_hiddenNodes(filenames, datafile_path, ofile):
    dist = distribution(filenames[0])
    if (dist=="uniform"):
        linestyle = '--'
    elif (dist=="normal"):
        linestyle = '-'
    for file in filenames:
        print "File:", file
        data = np.loadtxt(file, skiprows=1)
        cycles = data[1:, 0]
        energies = data[1:, 1]
        nv, nh = nNodes(file)
        plt.plot(cycles, energies, linestyle=linestyle, label="%d" %Decimal(nh))
    plt.xlabel("cycles")
    plt.ylabel(r"$E_L$")
    plt.title("Energy dependece on learning rates")
    plt.legend()
    # Show the minor grid lines with very faint and almost transparent grey lines
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
    plt.tight_layout()
    plt.savefig(datafile_path+ofile+"_"+dist+'_.pdf')
    plt.show()



"""
**********************************************************

Plotting of local energies as a function of variables
(learning rate and number of hidden nodes)

**********************************************************
"""
#******************************************************************
# Energies vs. variables (learning rate and number hidden nodes),
# with distribution as legend
#******************************************************************

# Function for plotting energies vs. learningrate, save wrt. number hidden nodes
def plotEnergy_learningrate(filenames, datafile_path, ofile):
    dist = distribution(filenames[0])
    nv, nh = nNodes(filenames[0])
    if (dist=="uniform"):
        linestyle = '--'
        color = 'r'
    elif (dist=="normal"):
        linestyle = '-'
        color = 'b'
    for file in filenames:
        print "File:", file
        data = np.loadtxt(file, skiprows=1)
        learningrates = data[1:, 0]
        energies = data[1:, 1]
        plt.semilogx(learningrates, energies, marker='+', color=color, markersize='7', label=dist, linestyle='none')
    plt.xlabel(r"$\eta$")
    plt.ylabel(r"$E_L$")
    plt.title(r"Relation between energy and learning rate, $\eta$")
    plt.legend()
    # Show the minor grid lines with very faint and almost transparent grey lines
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
    plt.tight_layout()
    plt.savefig(datafile_path+ofile+"_"+dist+'_nh_'+str(nh)+'_.pdf')
    plt.show()

# Function for plotting energies vs. number of hidden nodes, save wrt. learning rate
def plotEnergy_hiddenNodes(filenames, datafile_path, ofile):
    dist = distribution(filenames[0])
    learningrate = learningRate(filenames[0])
    if (dist=="uniform"):
        linestyle = '--'
        color = 'r'
    elif (dist=="normal"):
        linestyle = '-'
        color = 'b'
    for file in filenames:
        print "File:", file
        data = np.loadtxt(file, skiprows=1)
        nh = data[1:, 0]
        energies = data[1:, 1]
        plt.plot(nh, energies, marker='+', color=color, markersize='7', label=dist, linestyle='none')
    plt.xlabel(r"Number hidden nodes")
    plt.ylabel(r"$E_L$")
    plt.title("Energy as a function of hidden nodes")
    plt.legend()
    # Show the minor grid lines with very faint and almost transparent grey lines
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
    plt.tight_layout()
    plt.savefig(datafile_path+ofile+"_"+dist+'_eta_'+str(learningrate)+'_.pdf')
    plt.show()

"""
# Plot variance as a function of learningrate
filenames = sorted(glob.glob("../data/plots_and_data/b/RBM/uniform_vs_normal*.txt"))
print filenames
ofile = "variance_vs_learningrate.pdf"
learningRates = []
variances_uniform = []
variances_normal = []
for file in filenames:
    eta = float(file.split('_')[6])
    learningRates.append(eta)
    with open(file, 'r') as f:
        next(f)
        reader = csv.reader(f, delimiter='\t')
        for line in reader:
            if line[0]=='uniform':
                variances_uniform.append(float(line[2]))
            elif line[0]=='normal':
                variances_normal.append(float(line[2]))

plt.semilogx(learningRates, variances_uniform, 'bx', label='Uniform', markersize='5')
plt.semilogx(learningRates, variances_normal, 'rx', label='Normal', markersize='5')
plt.xlabel(r'$\eta$')
plt.ylabel(r'$\sigma^2$')
plt.title(r'Variance as a function of learning rate $\eta$')
plt.legend()
# Show the minor grid lines with very faint and almost transparent grey lines
#plt.minorticks_on()
#plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.grid('on')
plt.tight_layout()
plt.savefig(datafile_path+ofile)
plt.show()
"""
