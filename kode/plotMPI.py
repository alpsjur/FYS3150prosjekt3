import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

def plotVarT(ax, path, filename, nVar):
    data = np.loadtxt(path+filename)
    T = data[:,0]; var = data[:,nVar];
    #print(var)
    window_size = 5
    window = np.ones(int(window_size))/float(window_size)
    #var = np.convolve(var, window, "same")
    ax.plot(T[2:-2], var[2:-2], marker='.')

def findTc(path, filename, nVar):
    data = np.loadtxt(path+filename)
    T = data[:,0]; var = data[:,nVar];
    indexTc = np.argmax(var)
    return T[indexTc]


path = "../data/"
figdir = "../figurer/"
sns.set()
sns.set_style("darkgrid")
sns.set_palette("Set2")
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


variables = [r"$\left\langle E\right\rangle$",r"$C_V$",r"$\left\langle M\right\rangle$",r"$\chi$", r"$\left\langle M\right\rangle$"]
legends = ["L=40","L=60", "L=80", "L=100"]
L = [40,60,80,100]
names=["energy.pdf","heat_capacity.pdf","magnetization.pdf","susceptibility.pdf"]
files = ["metropolis40.dat", "metropolis60.dat", "metropolis80.dat", "metropolis100.dat", "metropolis60_zoom.dat"]
files_zoom = ["metropolis60_zoom.dat", "metropolis80_zoom.dat"]#["metropolis40_zoom.dat", "metropolis60_zoom.dat","metropolis80_zoom.dat", "metropolis100_zoom.dat"]

'''
for i in [1,2,3,4]:
    fig, ax = plt.subplots()
    for file in files:
        plotVarT(ax, path, file, i)
    ax.legend(legends,loc='upper center', bbox_to_anchor=(0.5, 1.18),fontsize=14, frameon=False, ncol=2)
    ax.set_xlabel("T",fontsize=14)
    ax.set_ylabel(variables[i-1],fontsize=14)
    #plt.savefig(figdir+names[i-1])
'''
for i in [1,2,3,4]:
    fig, ax = plt.subplots()
    for file in files_zoom:
        plotVarT(ax, path, file, i)
    ax.legend(legends,loc='upper center', bbox_to_anchor=(0.5, 1.18),fontsize=14, frameon=False, ncol=2)
    ax.set_xlabel("T",fontsize=14)
    ax.set_ylabel(variables[i-1],fontsize=14)


plt.show()

Tc = []
Tc_inf = []
for i in range(4):
    Tc.append(findTc(path, files[i], 2))
    print(findTc(path, files[i], 2))
    for j in range(i):
        Tc_inf.append((Tc[i]*L[i]-Tc[j]*L[j])/(L[i]-L[j]))
print(Tc_inf)
print(np.mean(Tc_inf))
