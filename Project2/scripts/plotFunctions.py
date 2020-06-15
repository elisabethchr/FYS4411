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
# filename = 2b-Initialization-x_RBMcycles_bruteForce_x_nv_x_nh_x_nMCsteps_x_nPart_x_nDim_x_stepL_x_eta_x_RBMcycles_x.txt
# where x represents an integer or float number

def initialization(filename):
    init = (filename.split('_')[0]).split('-')[2]
    return float(init)

def distribution(filename):
    dist = filename.split('_')[3]
    return dist

def solver(filename):
    solv = filename.split('_')[2]
    if solv=="bruteForce":
        solv = "Brute force"
    elif solv=="importance":
        solv = "Importance"
    elif solv=="gibbs":
        solv=="Gibbs"
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

#**************************************************************************
# Uniform vs normal weight distribution, with different initializations
#**************************************************************************

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
#    plt.title(r'Local energy, $\eta = %.1E$' %Decimal(eta))
    plt.legend()
    # Show the minor grid lines with very faint and almost transparent grey lines
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
    plt.tight_layout()
    plt.savefig(datafile_path+ofile+'_.pdf')
    plt.show()

    # header = "Distribution  <E_L>   sigma2  error"
    #
    # # write to file
    # dataset = np.array([distributions, mean, variance, error])
    # dataset = dataset.T
    # print dataset
    #
    # datafile = datafile_path+ofile+'_.txt'
    # with open(datafile, "w+") as datafile:
    #     np.savetxt(datafile, dataset, fmt = ['%s', '%s', '%s', '%s'], header=header, delimiter = '\t')

    return dataset


#*************************************************************
# Energies vs. RBM cycles for different distributions and
# distinct learning rate as title, with initialization
# intervals as legends
# (uniform or normal distributions)
#*************************************************************
def plotRBM_distributionIntervals(filenames, datafile_path, ofile):
    for file in filenames:
        dist = distribution(file)
        solv = solver(file)
        eta = learningRate(file)
        init = initialization(file)
        if (dist=="uniform"):
            linestyle = '--'
            label = "(-%.3f, %.3f)" %(init, init)
        elif (dist=="normal"):
            linestyle = '-'
            label = r"$\sigma$ = %.3f" %init
        print "File: ", file
        data = np.loadtxt(file, skiprows=1)
        cycles = data[1:, 0]
        energies = data[1:, 1]
        plt.plot(cycles, energies, linestyle=linestyle, label=label)
    plt.xlabel("cycles")
    plt.ylabel(r'$E_L$')
    plt.title(r"Energy dependence on %s distribution initializations, $\eta = $%.1E" %(dist, Decimal(eta)))
    plt.legend()
    # Show the minor grid lines with very faint and almost transparent grey lines
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
    plt.tight_layout()
    plt.savefig(datafile_path+ofile+"_"+solv+"_"+dist+'_.pdf')
    plt.show()

#****************************************************
# Energies vs. RBM cycles, variables as legends
# (learning rate and number of hidden nodes)
#****************************************************

# Function for plotting energies vs. RBM cycles for distinct distributions, with learning rate as legend
def plotRBM_learningrate(filenames, datafile_path, ofile):
    dist = distribution(filenames[0])
    for file in filenames:
        if (dist=="uniform"):
            linestyle = '--'
        elif (dist=="normal"):
            linestyle = '-'
        print "File:", file
        data = np.loadtxt(file, skiprows=1)
        cycles = data[1:, 0]
        energies = data[1:, 1]
        eta = learningRate(file)
        plt.plot(cycles, energies, linestyle=linestyle, label="%.1E" %Decimal(eta))
    plt.xlabel("cycles")
    plt.ylabel(r'$E_L$')
    plt.title("Energy dependence on learning rates")
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
        # if nh==8:
        #     print cycles
        #     print energies
        plt.plot(cycles, energies, linestyle=linestyle, label="%d" %Decimal(nh))
    plt.xlabel("cycles")
    plt.ylabel(r"$E_L$")
    plt.title("Energy dependence on hidden nodes")
    plt.legend()
    # Show the minor grid lines with very faint and almost transparent grey lines
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
    plt.tight_layout()
    plt.savefig(datafile_path+ofile+"_"+dist+'_.pdf')
    plt.show()


"""
**********************************************************

Plotting of local energies as a function of MC steps

**********************************************************
"""
#*********************************************************
# Energies vs. MC steps, with type of solver as legend
# (brute force, importance)
#*********************************************************

# Function for plotting energies vs. MC steps, with solvers as legends.
# Save wrt. learningrate, distribution, number hidden nodes and visible nodes,
# and number RBM cycles. PLotting every 1000'th step
# (MC steps are from the final optimization cycle)
def plotEnergy_MCsteps(filenames, datafile_path, ofile):
    RBM = RBMCycles(filenames[0])
    eta = learningRate(filenames[0])
    for file in filenames:
        print "File: ", file
        dist = distribution(file)
        nv, nh = nNodes(file)
        init = initialization(file)
        solv = solver(file)
        if (dist=="bruteForce"):
            linestyle = '--'
        elif (dist=="importance"):
            linestyle = '-'
        data = np.loadtxt(file, skiprows=1)
        step = data[:, 0]
        steps = step[0::1000]
        energy = data[:, 1]
        energies = energy[0::1000]
        plt.plot(steps, energies, linestyle='none', marker='.', label=solv)
    plt.legend()
    plt.xlabel("Metropolis steps")
    plt.ylabel(r"$E_L$")
    plt.title("Relation between energy and solver")
    plt.legend()
    plt.grid('on')
    plt.tight_layout()
    plt.savefig(datafile_path+ofile+"_"+dist+"_eta_"+str(eta)+"_nh_"+str(nh)+"_nv_"+str(nv)+"_RBM_"+str(RBM)+"_.pdf")
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

# Function for plotting energies vs. learningrate, with initializations as legends.
# Save wrt. number hidden nodes and number RBM cycles
def plotEnergy_learningrate(filenames, datafile_path, ofile):
    solv = solver(filenames[0])
    RBM = RBMCycles(filenames[0])
    for file in filenames:
        print "File:", file
        dist = distribution(file)
        nv, nh = nNodes(file)
        init = initialization(file)
        if (dist=="uniform"):
            linestyle = '--'
            label = "(-%.3f, %.3f)" %(init, init)
        elif (dist=="normal"):
            linestyle = '-'
            label = r"$\sigma=%.3f$" %init
        data = np.loadtxt(file, skiprows=1)
        learningrates = data[1:, 0]
        energies = data[1:, 1]
        plt.semilogx(learningrates, energies, marker='+', markersize='7', label=label, linestyle='none')
    plt.xlabel(r"$\eta$")
    plt.ylabel(r"$E_L$")
    plt.title(r"Relation between energy and learning rate, $\eta$ (%s distribution)" %dist)
    plt.legend()
    # Show the minor grid lines with very faint and almost transparent grey lines
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
    plt.tight_layout()
    plt.savefig(datafile_path+ofile+"_"+solv+"_"+dist+'_nh_'+str(nh)+"_RBM_"+str(RBM)+'_.pdf')
    plt.show()

# Function for plotting energies vs. number of hidden nodes, with initializations as legends
# Save wrt. learning rate
def plotEnergy_hiddenNodes(filenames, datafile_path, ofile):
    dist = distribution(filenames[0])
    learningrate = learningRate(filenames[0])
    RBM = RBMCycles(filenames[0])
    for file in filenames:
        init = initialization(file)
        if (dist=="uniform"):
            linestyle = '--'
            label = "(-%.3f, %.3f)" %(init, init)
        elif (dist=="normal"):
            linestyle = '-'
            label = r"$\sigma=%.3f$" %init
        print "File:", file
        data = np.loadtxt(file, skiprows=1)
        nh = data[1:, 0]
        energies = data[1:, 1]
        plt.plot(nh, energies, marker='+', markersize='7', label=label, linestyle='none')
    plt.xlabel(r"Number hidden nodes")
    plt.ylabel(r"$E_L$")
    plt.title("Energy as a function of hidden nodes")
    plt.legend()
    # Show the minor grid lines with very faint and almost transparent grey lines
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
    plt.tight_layout()
    plt.savefig(datafile_path+ofile+"_"+dist+'_eta_'+str(learningrate)+"_RBM_"+str(RBM)+'_.pdf')
    plt.show()


# Function for plotting energies vs. learning rates, with solvers as legends
# Save wrt. learning rate
def plotEnergy_learningrate_legendSolvers(filenames, datafile_path, ofile):
    dist = distribution(filenames[0])
    RBM = RBMCycles(filenames[0])
    nv, nh = nNodes(filenames[0])
    for file in filenames:
        global markerstyle
        markerstyle = '+'
        solv = solver(file)
        if (solv=="bruteForce"):
            markerstyle = '+'
        elif (solv=="Importance"):
            markerstyle = 'o'
        elif (solv=="Gibbs"):
            markerstyle = '.'
        data = np.loadtxt(file, skiprows=1)
        nh = data[1:, 0]
        energies= data[1:, 1]
        plt.plot(nh, energies, marker=markerstyle, linestyle='none', label=solv)
    plt.xlabel(r"$\eta$")
    plt.ylabel(r"$E_L$")
    plt.title("Energy as a function of learning rates")
    plt.legend()
    # Show the minor grid lines with very faint and almost transparent grey lines
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
    plt.tight_layout()
    plt.savefig(datafile_path+ofile+"_"+solv+"_"+dist+'_nh_'+str(nh)+"_RBM_"+str(RBM)+'_.pdf')
    plt.show()


# Function for plotting energies vs. number of hidden nodes, with solvers as legends
# Save wrt. learning rate
def plotEnergy_hiddenNodes_legendSolvers(filenames, datafile_path, ofile):
    dist = distribution(filenames[0])
    learningrate = learningRate(filenames[0])
    RBM = RBMCycles(filenames[0])
    for file in filenames:
        global markerstyle
        markerstyle = '+'
        solv = solver(file)
        if (solv=="bruteForce"):
            markerstyle = '+'
        elif (solv=="Importance"):
            markerstyle = 'o'
        elif (solv=="Gibbs"):
            markerstyle = '.'
        data = np.loadtxt(file, skiprows=1)
        nh = data[1:, 0]
        energies= data[1:, 1]
        plt.plot(nh, energies, marker=markerstyle, linestyle='none', label=solv)
    plt.xlabel(r"Number hidden nodes")
    plt.ylabel(r"$E_L$")
    plt.title("Energy as a function of hidden nodes")
    plt.legend()
    # Show the minor grid lines with very faint and almost transparent grey lines
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
    plt.tight_layout()
    plt.savefig(datafile_path+ofile+"_"+dist+'_eta_'+str(learningrate)+"_RBM_"+str(RBM)+'_.pdf')
    plt.show()
