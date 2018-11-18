import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import seaborn as sns

def plotRunningMean(ax, path, filename, nVar):
    data = np.loadtxt(path+filename)
    T = data[:,0]; var = data[:,nVar];
    #print(var)
    window_size = 11
    index = window_size//2
    window = np.ones(int(window_size))/float(window_size)
    var_con = np.convolve(var, window, "same")
    line, = ax.plot(T[index:-index], var_con[index:-index])
    linecolor=line.get_color()
    ax.plot(T,var, '.',color=linecolor)
    return linecolor

def plotVarT(ax, path, filename, nVar):
    data = np.loadtxt(path+filename)
    T = data[:,0]; var = data[:,nVar];
    #print(var)
    ax.plot(T,var, marker='.')

def findTc(path, filename, nVar):
    data = np.loadtxt(path+filename)
    T = data[:,0]; var = data[:,nVar];
    indexTc = np.argmax(var)
    return T[indexTc], var[indexTc]

def findTcInf(path, files, runningMean = False):
    Tc = []
    Tc_inf = []
    for i in range(len(files)):
        data = np.loadtxt(path+files[i])
        T = data[:,0]; var = data[:,2];
        if runningMean:
            window_size = 11
            window = np.ones(int(window_size))/float(window_size)
            var = np.convolve(var, window, "same")
        indexTc = np.argmax(var)
        Tc.append(T[indexTc])
        for j in range(i):
            print("i: ", i, "j: ", j, "Tc: ", (Tc[i]*L[i]-Tc[j]*L[j])/(L[i]-L[j]))
            Tc_inf.append((Tc[i]*L[i]-Tc[j]*L[j])/(L[i]-L[j]))
    return np.mean(Tc_inf), np.std(Tc_inf)


path = "../data/"
figdir = "../figurer/"
sns.set()
sns.set_style("darkgrid")
sns.set_palette("Set2")
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


variables = [r"$\left\langle E\right\rangle$",r"$C_V$",r"$\left\langle |M|\right\rangle$",r"$\chi$", r"$\left\langle M\right\rangle$"]
legends = ["L=40","L=60", "L=80", "L=100"]
L = [40,60,80,100]
names=["energy.pdf","heat_capacity.pdf","magnetization.pdf","susceptibility.pdf"]
files = ["metropolis40.dat", "metropolis60.dat", "metropolis80.dat", "metropolis100.dat"]
files_zoom = ["metropolis40_zoom.dat","metropolis60_zoom.dat", "metropolis80_zoom.dat", "metropolis100_zoom.dat"]

for i in [1,2,3,4]:
    fig, ax = plt.subplots()
    for file in files:
        plotVarT(ax, path, file, i)
    ax.legend(legends,loc='upper center', bbox_to_anchor=(0.5, 1.18),fontsize=14, frameon=False, ncol=2)
    ax.set_xlabel("T",fontsize=14)
    ax.set_ylabel(variables[i-1],fontsize=14)
    #plt.savefig(figdir+names[i-1])


'''
for i in [2]:
    fig, ax = plt.subplots()
    colors = []
    for file in files_zoom:
        colors.append(plotRunningMean(ax, path, file, i))
    #ax.legend(legends,loc='upper center', bbox_to_anchor=(0.5, 1.18),fontsize=14, frameon=False, ncol=2)
    ax.set_xlabel("T",fontsize=14)
    ax.set_ylabel(variables[i-1],fontsize=14)
    custom_lines = [Line2D([0], [0], color=colors[0]),
                    Line2D([0], [0], color=colors[1]),
                    Line2D([0], [0], color=colors[2]),
                    Line2D([0], [0], color=colors[3])]
    ax.legend(custom_lines,legends,loc='upper center', bbox_to_anchor=(0.5, 1.18),fontsize=14, frameon=False, ncol=2)
    #plt.savefig(figdir+"zoom.pdf")


Tc, std = findTcInf(path, files)
print("Forventningsverdi Tc: ", Tc, "Standardavvik: ", std)

Tc_zoom, std_zoom = findTcInf(path, files_zoom)
print("Forventningsverdi Tc: ", Tc_zoom, "Standardavvik: ", std_zoom)

Tc_zoom, std_zoom = findTcInf(path, files_zoom, runningMean = True)
print("Forventningsverdi Tc: ", Tc_zoom, "Standardavvik: ", std_zoom)

'''


plt.show()
