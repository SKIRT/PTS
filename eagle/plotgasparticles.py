#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.eagle.plotgasparticles Plot the gas particles for an EAGLE SKIRT-run in 3D.

# The facilities in this module serve to create a 3D plot of the gas particles for
# a particular EAGLE SKIRT-run, indicating star forming, cold and hot gas in different
# colors. The data is obtained from the SPH particles in the gas input file.

# ----------------------------------------------------------------------

# use a non-interactive back-end to generate high-quality vector graphics
import matplotlib
import matplotlib.pyplot as plt

# import 3D capabilities
from mpl_toolkits.mplot3d import Axes3D

# import standard modules
import os.path
import numpy as np

# import pts modules
from ..core.tools import archive as arch

# ----------------------------------------------------------------------

# define global parameters
extent = 30     # 3D box half-size in kpc
Tmax = 8000     # temperature cutoff for cold gas in K

# define corner points (we need to plot these to fix the extent of plotted volume)
cornerx = ( -extent, -extent, -extent, -extent,  extent,  extent,  extent,  extent )
cornery = ( -extent, -extent,  extent,  extent, -extent, -extent,  extent,  extent )
cornerz = ( -extent,  extent, -extent,  extent, -extent,  extent, -extent,  extent )

# ----------------------------------------------------------------------

# load columns text file in given directory and with name ending with given extension
def loadfile(inpath, extension):
    filenames = arch.listdir(inpath, extension)
    if len(filenames)!=1: raise ValueError("input file not found")
    filepath = os.path.join(inpath, filenames[0])
    return np.loadtxt(arch.opentext(filepath), unpack=True)

# ----------------------------------------------------------------------

## This function creates a 3D plot of the gas particles for  a particular EAGLE SKIRT-run,
# indicating star forming, cold and hot gas in different colors.
# The data is obtained from the SPH particles in the gas input file.
# The output plot is placed in the SKIRT-run's visualization directory.
def plotgasparticles(skirtrun):
    # load the data
    x,y,z,h,M,Z,T = loadfile(skirtrun.inpath(), "_gas.dat")

    # convert positions from pc to kpc
    x /= 1000.
    y /= 1000.
    z /= 1000.

    # create masks for different particle sets
    inbox = (np.abs(x)<extent) & (np.abs(y)<extent) & (np.abs(z)<extent)
    sf = T<=0
    cold = (0<T) & (T<=Tmax)
    hot = T>Tmax

    # setup the 3D figure
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_aspect('equal')
    ax.set_xlim(-extent,extent)
    ax.set_ylim(-extent,extent)
    ax.set_zlim(-extent,extent)
    ax.view_init(elev=20, azim=45)
    ax.set_xlabel("x (kpc)", fontsize='medium')
    ax.set_ylabel("y (kpc)", fontsize='medium')
    ax.set_zlabel("z (kpc)", fontsize='medium')
    ax.set_title("runid {} -- {}".format(skirtrun.runid(), skirtrun.prefix()), fontsize='medium')

    # add the scatter plots
    ax.scatter(x[inbox&sf], y[inbox&sf], z[inbox&sf], marker=".", lw=0, color='b')
    ax.scatter(x[inbox&cold], y[inbox&cold], z[inbox&cold], marker=".", lw=0, color='c')
    ax.scatter(x[inbox&hot], y[inbox&hot], z[inbox&hot], marker=".", lw=0, color='r')
    ax.scatter(cornerx, cornery, cornerz, marker=".", lw=0, color='w')

    # add the legend (through a work-around)
    proxy1 = matplotlib.lines.Line2D([0],[0], linestyle='none', c='b', marker = '.')
    proxy2 = matplotlib.lines.Line2D([0],[0], linestyle='none', c='c', marker = '.')
    proxy3 = matplotlib.lines.Line2D([0],[0], linestyle='none', c='r', marker = '.')
    ax.legend([proxy1, proxy2, proxy3], ["star forming gas", "cold gas", "hot gas"], numpoints = 1,
                loc='upper right', prop={'size':'small'})

    # save the figure
    plotpath = os.path.join(skirtrun.vispath(), skirtrun.prefix()+"_gas_particles.pdf")
    plt.savefig(plotpath, bbox_inches='tight', pad_inches=0.25)
    plt.close()
    print "Created PDF plot file " + plotpath

# ----------------------------------------------------------------------
