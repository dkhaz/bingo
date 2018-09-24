######################################################################
# Python code to produce 2D plots
# Written by V. Sreenath and L. Sriramkumar 
######################################################################
# Importing libraries
import matplotlib
from numpy import *
from pylab import *
import pylab
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from matplotlib import rc, rcParams, pyplot
from matplotlib.pyplot import figure, axes, plot, xlabel, ylabel, title, grid, savefig, show
from matplotlib.ticker import MaxNLocator
from pyx import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
######################################################################
figure(1,figsize=(8,6))
######################################################################
rc('font', family='TimesRoman', weight = 'extra bold', size = 18.0)
rc('text', usetex=True)
rc('axes', linewidth = 2, labelsize = 'large')  
rc('xtick', labelsize= 'medium')
rcParams['xtick.major.size'] = 8.0 
rcParams['xtick.minor.size'] = 4.0
rcParams['xtick.major.pad'] = 8.0 
rcParams['xtick.minor.pad'] = 8.0
rc('ytick', labelsize= 'medium')  
rcParams['ytick.major.size'] = 8.0 
rcParams['ytick.minor.size'] = 0.0
rcParams['ytick.major.pad'] = 8.0 
rcParams['ytick.minor.pad'] = 8.0
rc('lines', linewidth = 1.5, markeredgewidth=1.5)
rc('savefig', dpi=300)
######################################################################
# Fixing space around the plot
pylab.axes([0.125,0.125,0.825,0.825])
######################################################################
#Reading data from files and plotting them
#read from textfile to x and y
data = genfromtxt('plots/F_nl_2d.txt')
x = data[:,0]
y = data[:,1]
z = data[:,2]
xi = np.linspace(0,1,400)
yi = np.linspace(0.5,1,400)
zi = griddata(x,y,z,xi,yi,interp='linear')
f = pyplot.figure()
ax = f.gca()
CS = plt.contourf(xi,yi,zi,50,cmap=plt.cm.jet)
plt.axis([0, 1, 0.5, 1])
pylab.xlabel(r'$k_3/k_1$')
pylab.ylabel(r'$k_2/k_1$')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right",size="3%", pad=0.25)
m1=round(zi.min(),3)
m2=round(zi.max(),3)
md=round(m1+(m2-m1)/2,3)
print m1,m2,md
cbar=plt.colorbar(ticks=[m1,md,m2], cax=cax)
cbar.set_ticklabels([m1,md,m2]) 
ax.set_aspect(1)
######################################################################
pylab.savefig('plots/F_nl_2d.eps',bbox_inches='tight')
#pylab.show()
###################################################################### 