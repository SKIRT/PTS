#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.scaleeagle Plot scaling relations for EAGLE galaxies.
#
# This script plots scaling relations for a set of EAGLE galaxies using SKIRT simulation results for those galaxies.
# For each galaxy, a panchromatic SKIRT simulation must have produced at least the following output files:
# "log.txt", "parameters.xml", "wavelengths.dat", "luminosities.dat", "ds_convergence.dat", and "instrname_sed.dat"
# for one or more instruments.
# The file path and other parameters are hardcoded in the script.

# -----------------------------------------------------------------

# the path of the directory containing the SKIRT simulation results for the galaxies to be processed
eagle_skirt_path = "/Users/pcamps/EAGLE/Snapshot50/sed"

# -----------------------------------------------------------------

# import standard modules
import os.path
import numpy as np

# use a non-interactive back-end to generate high-quality vector graphics
import matplotlib
if matplotlib.get_backend().lower() != "pdf": matplotlib.use("pdf")
import matplotlib.pyplot as plt

# import relevant pts modules
from pts.skirtsimulation import createsimulations
from pts.skirtsimulation import Msun
from pts.skirtsimulation import pc

# -----------------------------------------------------------------

# dictionary of tuples representing a type of plot axis
#  key: axis type identifier
#  value : ( label, lambda function for value given simulation and instrument name )
axisdefinitions = {

    # intrinsic properties
    "logMstar" : ( r"$\log_{10}(M_*\,[M_\odot])$", lambda sim, name: np.log10(sim.starmass()) ),
    "logMdust" : ( r"$\log_{10}(M_{dust}\,[M_\odot])$", lambda sim, name: np.log10(sim.dustmass()) ),
    "logMdust/Mstar" : ( r"$\log_{10}(M_{dust}/M_*)$", lambda sim, name: np.log10(sim.dustmass()/sim.starmass()) ),

    # colors
    "g-r" : ( r"$g-r$", lambda sim, name: sim.color(name,"g-r") ),
    "NUV-r" : ( r"$\mathrm{NUV}-r$", lambda sim, name: sim.color(name,"NUV-r") ),

    # magnitudes
    "r" : ( r"$M_r$", lambda sim, name: sim.absmagnitude(name,"r") ),

    # flux densities (in Jansky)
    "logS250" : ( r"$\log_{10}(S_{250}\,[\mathrm{Jy}])$", lambda sim, name: np.log10(sim.fluxdensityFreq(name,250)) ),
    "logS500" : ( r"$\log_{10}(S_{500}\,[\mathrm{Jy}])$", lambda sim, name: np.log10(sim.fluxdensityFreq(name,500)) ),

    # ratios of flux densities (in Jansky)
    "S250/S350" : ( r"$S_{250}/S_{350}$", lambda sim, name: sim.fluxdensityFreq(name,250) / sim.fluxdensityFreq(name,350) ),
    "S350/S500" : ( r"$S_{350}/S_{500}$", lambda sim, name: sim.fluxdensityFreq(name,350) / sim.fluxdensityFreq(name,500) ),
    "f250/f500" : ( r"$f_{250}/f_{500}$", lambda sim, name: sim.fluxdensityFreq(name,250) / sim.fluxdensityFreq(name,500) ),

    # luminosity densities (based on flux in Jansky)
    "logL250" : ( r"$\log_{10}(L_{250}\,[\mathrm{W/Hz}])$", lambda sim, name: np.log10(sim.luminositydensityFreq(name,250)) ),
    "logL500" : ( r"$\log_{10}(L_{500}\,[\mathrm{W/Hz}])$", lambda sim, name: np.log10(sim.luminositydensityFreq(name,500)) ),
    "logL250/LNUV" : ( r"$\log_{10}(L_{250}/L_\mathrm{NUV})$", lambda sim, name:
                        np.log10(sim.fluxdensityFreq(name,250) / sim.fluxdensityFreq(name,"NUV")) ),

    # integrated luminosities
    "logLdust" : ( r"$\log_{10}(L_{dust}\,[L_\odot])$", lambda sim, name:
                        np.log10(sim.luminosityforflux(sim.integratedflux(name,8,1000))) ),

    # other ratios
    "logMdust/f350" : ( r"$\log_{10}(M_{dust}/f_{350}/D^2\,[\mathrm{kg}\,\mathrm{W}^{-1}\,\mathrm{Hz}])$", lambda sim, name:
                        np.log10( (sim.dustmass()*Msun) /
                                  (sim.fluxdensityFreq(name,350)*1e-26) / (sim.instrumentdistance()*pc)**2 ) ),
    "Mgas/Mtot" : ( r"$M_{gas}/(M_*+M_{gas})$", lambda sim, name: sim.gasmass() / (sim.starmass()+sim.gasmass()) ),

}

# -----------------------------------------------------------------

# this function produces a one-page pdf file with one or more scaling plots for a set of SKIRT-EAGLE simulations.
# arguments:
# - simdir: relative or absolute path to the directory containing the output of the SKIRT simulations to be plotted
# - plotfile: name of the output file including the .pdf extension; may include a relative or absolute path
# - plotdefs: sequence of plot definitions; each item is a 2-tuple of axis type identifiers (key in above dict)
# - pagesize: a 2-tuple specifying the size of the complete page in inch
# - layout: a 2-tuple specifying the number of columns and rows in the layout of the plots
#           (the layout must accomodate all items in the plotdefs sequence)
def plotscaling(simdir, plotfile, plotdefs, layout=(2,3), pagesize=(8,12)):
    simdir = os.path.realpath(os.path.expanduser(simdir))
    plotfile = os.path.realpath(os.path.expanduser(plotfile))
    assert plotfile.endswith(".pdf")

    # create a simulation object for each log file in the directory
    simulations = createsimulations(eagle_skirt_path)
    print "Plotting", len(plotdefs), "scaling relations for", len(simulations), "galaxies"

    # setup the figure
    figure = plt.figure(figsize=pagesize)
    figure.subplots_adjust(wspace=0.35, hspace=0.25)

    # loop over the plots
    plotindex = 0
    for xaxis,yaxis in plotdefs:
        xlabel,xvalue = axisdefinitions[xaxis]
        ylabel,yvalue = axisdefinitions[yaxis]

        # start the appropriate subplot
        plotindex += 1
        plt.subplot(layout[1], layout[0], plotindex)

        # setup the x and y values for each of the simulations
        x = [ ]
        y = [ ]
        for sim in simulations:
            for name in sim.instrumentnames():
                x += [ xvalue(sim,name) ]
                y += [ yvalue(sim,name) ]

        # plot the scaling relation
        plt.scatter(x, y)

        # add axis labels
        plt.xlabel(xlabel, fontsize='large')
        plt.ylabel(ylabel, fontsize='large')

    # save and close the figure
    plt.savefig(plotfile, bbox_inches='tight', pad_inches=0.25)
    plt.close()
    print "Created plot file", plotfile

# -----------------------------------------------------------------

# create some specific plots

plotscaling(eagle_skirt_path, "scaling1.pdf",
            [ ("r","g-r"), ("r","NUV-r"), ("r","logS250"), ("r","logS500"), ("r","logL250"), ("r","logL500") ] )

plotscaling(eagle_skirt_path, "scaling2.pdf",
            [ ("S250/S350","S350/S500"), ("r","r"),
              ("logMstar","logL250"), ("logMstar","logL500"),
              ("logMstar","logMdust"), ("logMstar","logMdust/Mstar") ] )

plotscaling(eagle_skirt_path, "scaling3.pdf",
            [ ("logMstar","logLdust"), ("logMstar","logL250/LNUV"),
              ("r","r"), ("r","r"),
              ("r","r"), ("r","r") ] )

plotscaling(eagle_skirt_path, "scaling4.pdf",
            [ ("logMdust/f350","f250/f500"), ("logMstar","logMdust/Mstar"),
              ("NUV-r","logMdust/Mstar"), ("Mgas/Mtot","logMdust/Mstar"),
              ("r","r"), ("r","r") ] )

# -----------------------------------------------------------------
