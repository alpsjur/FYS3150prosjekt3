'''
Kode for å plotte fra dataen generert av koden som ikke er paralellisert
'''

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.mlab as mlab
import numpy as np
import seaborn as sns


def plotDatavsMCS(ax, path, filename, i, line='.'):
    data = np.fromfile(path+filename)
    mcs = data[0::8]; y = data[i::8]
    #ax.set_xscale('log')
    line, = ax.plot(mcs, y, line, marker='.')
    return line

def plotProbability(ax, path, filename):
    data = np.loadtxt(path+filename)
    #henter ut forventningverdi og standardavvik
    expE = data[0,0]; var = data[0,1]; sigma = np.sqrt(var/20/20)
    #heter ut energiene og tilhørende sannsynlighet
    E = data[1:,0]; p = data[1:,1]

    #følgende kode er for å plotte standardavvik, forventningsverdi og normalfordeling
    #ax.axvspan(expE-sigma, expE+sigma, alpha=0.4,color='lightblue') #standard
    #ax.axvline(x=expE,color="grey")
    #gaus = mlab.normpdf(E, expE, sigma)
    #ax.plot(E,gaus/np.sum(gaus))

    line, = ax.plot(E,p,'.')
    #ax.bar(E, p,0.01)
    return line



path = "../data/"
figdir = "../figurer/"
sns.set(color_codes=True)
sns.set_style("darkgrid")
sns.set_palette("Set2")
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


#plotter energi og magnetisering mot mcs
fig, ax = plt.subplots(2,1, sharex = True)
line1 = plotDatavsMCS(ax[0], path, "EMvsMCS_unordered_1.bin", 2)
line2 = plotDatavsMCS(ax[0], path, "EMvsMCS_ordered_1.bin", 2)

line3 = plotDatavsMCS(ax[0], path, "EMvsMCS_unordered_24.bin", 2)
line4 = plotDatavsMCS(ax[0], path, "EMvsMCS_ordered_24.bin", 2)

plotDatavsMCS(ax[1], path, "EMvsMCS_unordered_1.bin", 4)
plotDatavsMCS(ax[1], path, "EMvsMCS_ordered_1.bin", 4)

plotDatavsMCS(ax[1], path, "EMvsMCS_unordered_24.bin", 4)
plotDatavsMCS(ax[1], path, "EMvsMCS_ordered_24.bin", 4)

custom_lines = [Line2D([0], [0], color=line1.get_color()),
                Line2D([0], [0], color=line2.get_color()),
                Line2D([0], [0], color=line3.get_color()),
                Line2D([0], [0], color=line4.get_color())]

ax[0].legend(custom_lines,["Unordered T=1.0","Ordered T=1.0", "Unordered T=2.4", "Ordered T=2.4"],\
              loc='upper center', bbox_to_anchor=(0.5, 1.4),fontsize=14, frameon=False, ncol=2)
ax[1].set_xlabel("mcs",fontsize=14)
ax[0].set_ylabel(r"$\left\langle E\right\rangle$",fontsize=14)
ax[1].set_ylabel(r"$\left\langle |M|\right\rangle$",fontsize=14)
#plt.savefig(figdir+"EMvsMCS.pdf")

#plotter antall aksepterte nye konfigurasjoner mot mcs
fig, ax = plt.subplots()
plotDatavsMCS(ax, path, "count_unordered_1.bin", 7, '-')
plotDatavsMCS(ax, path, "count_ordered_1.bin", 7, '-')
plotDatavsMCS(ax, path, "count_unordered_24.bin", 7, '-')
plotDatavsMCS(ax, path, "count_ordered_24.bin", 7, '-')
#mcs = np.linspace(10,10**6)
#nSpins = 20*20
#ax.plot(mcs, nSpins*mcs, '--', color="lightgrey")

ax.set_xscale('log')
ax.set_yscale('log')

ax.legend(custom_lines,["Unordered T=1.0","Ordered T=1.0", "Unordered T=2.4", "Ordered T=2.4"],\
              loc='upper center', bbox_to_anchor=(0.5, 1.18),fontsize=14, frameon=False, ncol=2)

ax.set_xlabel("mcs",fontsize=14)
ax.set_ylabel("\# accepted configurations",fontsize=14)
#plt.savefig(figdir+"nAccepted.pdf")


#plotter sannsynlighetsfordelingen
fig, ax = plt.subplots(2, 1, sharex=True)
ax[0].plot([-2],[0])
line1 = plotProbability(ax[0], path, "probability1.dat")
line2 = plotProbability(ax[1], path, "probability24.dat")

ax[1].set_xlabel(r"E", fontsize=14)
#fig.text(0.5, 0.03, 'E',  ha='center',fontsize=14)
fig.text(0.02, 0.5, r'p(E)',  va='center', rotation='vertical',fontsize=14)

custom_lines = [Line2D([0], [0], color=line1.get_color()),
                Line2D([0], [0], color=line2.get_color())]

ax[0].legend(custom_lines,["T=1.0","T=2.4"], loc='upper center',
             bbox_to_anchor=(0.5, 1.3),fontsize=14, frameon=False, ncol=2)
#plt.savefig(figdir+"probability.pdf")

plt.show()
