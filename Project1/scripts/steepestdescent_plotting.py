import numpy as np
import matplotlib.pyplot as plt
import sys, glob, os

filenames = sorted(glob.glob('../data/f/*131072_.txt'), key=os.path.getmtime)

def alpha(filename):
    var = filename.split("_")
    return var[2]

def Steps(filename):
    var = filename.split("_")
    return var[6]

iter = []
alphas = []
E = []
legends = []
final_alphas = []
nSteps = []
for i in filenames:
    print i
    data = np.loadtxt(i, skiprows=1)
    iter = data[:, 0]
    alphas = data[:, 1]
    E = data[:, 2]
    a = alpha(i)
    legends.append(a)
    final_alphas.append(alphas[-1])
    nSteps.append(Steps(i))

    plt.plot(iter, alphas, linestyle="--")

plt.xlabel(r'Iterations')
plt.ylabel(r'$\alpha$')
plt.title("Steepest descent")
plt.legend(legends)
plt.grid(b=True, which='major', color='#666666', linestyle='-')
# Show the minor grid lines with very faint and almost transparent grey lines
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.tight_layout()
plt.grid('on')
#plt.savefig('../data/plots/f/steepest_descent_nSteps_'+nSteps[0]+'.pdf')
plt.show()

mean_alpha = np.mean(final_alphas)
std_alpha = np.std(final_alphas)
print mean_alpha, std_alpha
