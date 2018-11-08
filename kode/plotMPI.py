import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

def plotVarT(ax, path, filename, nVar):
    data = np.fromfile(path+filename)
    T = data[0::6]; var = data[nVar::6];
    ax.plot(T, var)

path = "../data/"
figdir = "../figurer/"
sns.set()
sns.set_style("darkgrid")
sns.set_palette("Set2")
plt.rc('text', usetex=True)
plt.rc('font', family='serif')




fig, ax = plt.subplots()
plotVarT(ax, path, "metropolis40.bin", 2)
plt.show()
