import matplotlib.pyplot as plt
from scipy.stats import norm
import numpy as np
import seaborn as sns


def plotDatavsMCS(ax, path, filename, i):
    data = np.loadtxt(path+filename)
    mcs = data[:,0]; y = data[:,i]
    ax.set_xscale('log')
    ax.plot(mcs, y, marker='.')

def plotProbability(ax, path, filename):
    data = np.loadtxt(path+filename)
    expE = data[0,0]; sigma = np.sqrt(data[0,1])
    E = data[1:,0]; p = data[1:,1]
    #ax.axvspan(expE-sigma, expE+sigma, alpha=0.4,color='lightblue')
    ax.axvline(x=expE,color="grey")
    ax.plot(E,p,'.')



path = "../data/"
figdir = "../figurer/"
sns.set()
sns.set_style("darkgrid")
sns.set_palette("Paired")
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

'''
fig, ax = plt.subplots(2,1, sharex = True)
plotDatavsMCS(ax[0], path, "EMvsMCS_ordered_1.dat", 2)
plotDatavsMCS(ax[0], path, "EMvsMCS_unordered_1.dat", 2)

plotDatavsMCS(ax[0], path, "EMvsMCS_ordered_24.dat", 2)
plotDatavsMCS(ax[0], path, "EMvsMCS_unordered_24.dat", 2)


plotDatavsMCS(ax[1], path, "EMvsMCS_ordered_1.dat", 4)
plotDatavsMCS(ax[1], path, "EMvsMCS_unordered_1.dat", 4)


plotDatavsMCS(ax[1], path, "EMvsMCS_ordered_24.dat", 4)
plotDatavsMCS(ax[1], path, "EMvsMCS_unordered_24.dat", 4)

ax[0].legend(["Ordered T=1.0","Unordered T=1.0", "Ordered T=2.4", "Unordered T=2.4"],\
              loc='upper center', bbox_to_anchor=(0.5, 1.4),fontsize=14, frameon=False, ncol=2)
ax[1].set_xlabel("mcs",fontsize=14)
ax[0].set_ylabel(r"$<E>$",fontsize=14)
ax[1].set_ylabel(r"$<|M|>$",fontsize=14)
#plt.savefig(figdir+"EMvsMCS.pdf")


fig, ax = plt.subplots()
plotDatavsMCS(ax, path, "EMvsMCS_ordered_1.dat", -1)
plotDatavsMCS(ax, path, "EMvsMCS_unordered_1.dat", -1)
plotDatavsMCS(ax, path, "EMvsMCS_ordered_24.dat", -1)
plotDatavsMCS(ax, path, "EMvsMCS_unordered_24.dat", -1)

ax.legend(["Ordered T=1.0","Unordered T=1.0", "Ordered T=2.4", "Unordered T=2.4"],\
              loc='upper center', bbox_to_anchor=(0.5, 1.18),fontsize=14, frameon=False, ncol=2)
ax.set_yscale('log')
ax.set_xlabel("mcs",fontsize=14)
ax.set_ylabel("\# accepted configurations",fontsize=14)
#plt.savefig(figdir+"nAccepted.pdf")
'''
sns.set_palette("Set2")

fig, ax = plt.subplots(2, 1, sharex=True)
plotProbability(ax[0], path, "probability1.dat")
plotProbability(ax[1], path, "probability24.dat")

plt.show()
