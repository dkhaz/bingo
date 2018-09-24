######################################################################
# Python code to produce 3D plots
# Written by V. Sreenath and L. Sriramkumar 
######################################################################
# CAUTION : Contour values need to be changed for different runs. Here 
# we provide a specific case of quadratic potential with step for 
# a particular best fit potential parameters.
######################################################################

import scipy as np
from pylab import *
from mayavi import mlab
#Read in data
data = np.genfromtxt('plots/F_nl_3d.txt')
a = data[:,0]
b = data[:,1]
c = data[:,2]
f = data[:,3]
print f.min(), f.max()

#Reshape data to match the format
# xyz data should be a dim x dim x dim matrix
dim = 60
a1 = a.reshape(dim,dim,dim)
b1 = b.reshape(dim,dim,dim)
c1 = c.reshape(dim,dim,dim)
values = f.reshape(dim,dim,dim)

# Specifies background and forground colour and size of figure
mlab.figure(1, size=(500, 250), fgcolor=(0, 0, 0), bgcolor=(1,1,1))

# Plots isosurfaces and specifies contours and opacity
mlab.contour3d(a1, b1, c1, values, contours=[-3.2926,-3.0, -2.0, -1.0,1.0, 2.0, 3.0, 4.1413], transparent=False,opacity=0.9, colormap='jet')

#Colorbar and its properties
cb=mlab.colorbar(title='', orientation='vertical', nb_labels=10)
cb.label_text_property.font_family = 'times'
cb.label_text_property.font_size = 10

#Axes and its properties
mlab.outline()
ax=mlab.axes(xlabel= '', ylabel='', zlabel='', ranges=(1.0e-5,1.0e-2,1.0e-5,1.0e-2,1.0e-5,1.0e-2))
ax.axes.font_factor=0.5
ax.label_text_property.font_family = 'times'

mlab.show()

 