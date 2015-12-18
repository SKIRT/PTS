#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.eagle.plotdensitycurves Plot stellar, gas and dust
# densities in function of the galaxy radius for an EAGLE SKIRT-run.
#
# The facilities in this module serve to plot stellar, gas and dust
# densities in function of the galaxy radius for a particular EAGLE SKIRT-run.
# The data is calculated based on the SPH particles in the input files,
# assuming al mass is concentrated in a particle's center position.

# ----------------------------------------------------------------------

# use a non-interactive back-end to generate high-quality vector graphics
import matplotlib.pyplot as plt

# import standard modules
import os.path
import numpy as np

# import pts modules
from ..core.tools import archive as arch

# ----------------------------------------------------------------------

# load columns text file in given directory and with name ending with given extension
def loadfile(inpath, extension):
    filenames = arch.listdir(inpath, extension)
    if len(filenames)!=1: raise ValueError("input file not found")
    filepath = os.path.join(inpath, filenames[0])
    return np.loadtxt(arch.opentext(filepath), unpack=True)

# ----------------------------------------------------------------------

## This function creates a PDF plot with histograms of stellar, gas and dust
# densities in function of the galaxy radius for a particular EAGLE SKIRT-run.
# The data is calculated based on the SPH particles in the input files,
# assuming al mass is concentrated in a particle's center position.
# The output plot is placed in the SKIRT-run's visualization directory.
def plotdensitycurves(skirtrun):
    # setup the figure
    figure = plt.figure(figsize=(10,6))
    rmax = 50  # kpc

    # load and plot the stars
    x,y,z,h,M,Z,t = loadfile(skirtrun.inpath(), "_stars.dat")
    r = np.sqrt(x*x + y*y + z*z)/1000  # kpc
    r[r>rmax] = rmax
    plt.hist(r, weights=M, bins=25, range=(0,rmax), histtype='step', log=True, color='b', label="stars")

    # load and plot the gas
    x,y,z,h,Mgas,Z,T = loadfile(skirtrun.inpath(), "_gas.dat")
    r = np.sqrt(x*x + y*y + z*z)/1000  # kpc
    r[r>rmax] = rmax
    M = Mgas.copy()
    M[np.abs(T)>75000] = 0
    if np.any(M):
        plt.hist(r, weights=M*Z, bins=25, range=(0,rmax), histtype='step', log=True, color='m', ls='dashed', label="metals (T<75000K)")
    M = Mgas.copy()
    M[np.abs(T)>8000] = 0
    if np.any(M):
        plt.hist(r, weights=M*Z, bins=25, range=(0,rmax), histtype='step', log=True, color='m', ls='dotted', label="metals (T<8000K)")
    M = Mgas.copy()
    M[T>8000] = 0
    if np.any(M):
        plt.hist(r, weights=M*Z, bins=25, range=(0,rmax), histtype='step', log=True, color='m', ls='solid', label="metals (T<8000K or SFR>0)")

    # load and plot the hii regions
    try:
        x,y,z,h,SFR,Z,logC,P,fPDR = loadfile(skirtrun.inpath(), "_hii.dat")
        r = np.sqrt(x*x + y*y + z*z)/1000  # kpc
        r[r>rmax] = rmax
        plt.hist(r, weights=SFR*1e7, bins=25, range=(0,rmax), histtype='step', log=True, color='c', label="hii regions")
    except ValueError:
        pass

    # add axis labels, legend and title
    plt.grid('on')
    plt.xlabel("r (kpc)", fontsize='medium')
    plt.ylabel("Mass (Msun)", fontsize='medium')
    plt.ylim(1e4, 1e9)
    plt.legend(loc='upper right', prop={'size':'small'})
    plt.title("runid {} -- {}".format(skirtrun.runid(), skirtrun.prefix()), fontsize='medium')

    # save the figure
    plotpath = os.path.join(skirtrun.vispath(), skirtrun.prefix()+"_density_curves.pdf")
    plt.savefig(plotpath, bbox_inches='tight', pad_inches=0.25)
    plt.close()
    print "Created PDF plot file " + plotpath

# ----------------------------------------------------------------------
