#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package eagle.plotgreybodyfit Plot a modified blackbody fit to
# the dust continuum emission for an EAGLE SKIRT-run.
#
# The facilities in this module serve to plot a modified blackbody fit to
# the dust continuum emission for a particular EAGLE SKIRT-run.

# ----------------------------------------------------------------------

# use a non-interactive back-end to generate high-quality vector graphics
import matplotlib.pyplot as plt

# import standard modules
import os.path
import numpy as np
from scipy.optimize import curve_fit

# import pts modules
from ..core.tools import archive as arch
from ..core.basics.filter import Filter

# ----------------------------------------------------------------------

# universal constants and units
c = 2.99792458e8        # m/s
h = 6.62606957e-34      # J s
k = 1.3806488e-23       # J/K
pc = 3.08567758e16      # m
Msun = 1.9891e30        # kg

# global configuration values
beta = 2
kappa350 = 0.192        # m2/kg (Cortese)
nu350 = c/350e-6        # Hz
D = 1e6 * pc            # m  (instrument distance in ski file)

# returns grey body flux in Jy for wavelength(s) specified in micron,
# for given temperature T (K) and dust mass M (Msun),
# with global beta, kappa, and distance
def greybody(wave, T, M):
    nu = c / (wave * 1e-6)                                  # Hz
    kappa = kappa350 * (nu/nu350)**beta                     # m2/kg
    Bnu = 2*h*nu**3/ c**2 / (np.exp((h*nu)/(k*T)) - 1)      # W/m2/Hz
    flux = M * Msun * kappa * Bnu / D**2                    # W/m2/Hz
    return flux * 1e26                                      # Jy

# ----------------------------------------------------------------------

# returns temperature T (K) and dust mass M (Msun) for best fit with given data points and uncertainties
def fitgreybody(lambdav, fluxv, sigmav):
    # optimize the fit
    popt, pcov = curve_fit(greybody, lambdav, fluxv, p0=(17,1e7), sigma=sigmav, absolute_sigma=False, maxfev=5000)
    return popt

# ----------------------------------------------------------------------

## This function creates a PDF plot of a modified blackbody fit to
# the dust continuum emission for a particular EAGLE SKIRT-run,
# also listing the corresponding temperature.
# The output plot is placed in the SKIRT-run's visualization directory.
def plotgreybodyfit(skirtrun):
    simulation = skirtrun.simulation()

    # setup the figure
    figure = plt.figure(figsize=(10,6))
    plt.xscale('log')
    plt.yscale('log')

    # load and plot the SED
    filepath = simulation.seddatpaths()[0]
    lambdav, fluxv = np.loadtxt(arch.opentext(filepath), usecols=(0,1), unpack=True)
    plot = (lambdav>=10) & (lambdav<=1000)
    plt.plot(lambdav[plot], fluxv[plot], color='b', label="SKIRT galaxy SED")

    # load and plot the Herschel continuum data points (160, 250, 350, 500 micron)
    info = { }
    infofile = arch.listdir(skirtrun.vispath(), "_info.txt")[0]
    for line in arch.opentext(os.path.join(skirtrun.vispath(),infofile)):
        if not line.startswith("#"):
            key,dummy,value = line.split(None, 2)
            info[key] = float(value)
    waves = np.array( [ Filter(fs).pivotwavelength() for fs in ("Pacs.red","SPIRE.PSW","SPIRE.PMW","SPIRE.PLW")] )
    fluxes = np.array(( info['instr_xy_fluxdensity_pacs_red_continuum'],
                        info['instr_xy_fluxdensity_spire_psw_continuum'],
                        info['instr_xy_fluxdensity_spire_pmw_continuum'],
                        info['instr_xy_fluxdensity_spire_plw_continuum']  ))
    sigmas = np.array(( 3,1,1,3 ))      # pacs is less sensitive; longer wavelength fluxes are harder to measure
    plt.scatter(waves, fluxes, color='r', marker='*', label="Mock PACS/SPIRE fluxes")

    # fit a grey body to the Herschel fluxes and plot the result
    T,M = fitgreybody(waves, fluxes, sigmas)
    plt.plot(lambdav[plot], greybody(lambdav[plot], T,M), color='m',
        label=r"Grey body fit $T={:.2f},\,M_\mathrm{{dust}}={:.2e}\,M_\odot$".format(T,M))

    # add axis labels, legend and title
    plt.grid('on')
    plt.xlabel(r"$\lambda\,(\mu \mathrm{m})$", fontsize='medium')
    plt.ylabel(simulation.fluxlabel(), fontsize='medium')
    plt.xlim(10, 1000)
    ymax = fluxv[plot].max()
    plt.ylim(ymax*1.1e-3, ymax*1.1)
    plt.legend(loc='upper left', prop={'size':'small'})
    plt.title("runid {} -- {}".format(skirtrun.runid(), skirtrun.prefix()), fontsize='medium')

    # save the figure
    plotpath = os.path.join(skirtrun.vispath(), skirtrun.prefix()+"_dust_body_fit.pdf")
    plt.savefig(plotpath, bbox_inches='tight', pad_inches=0.25)
    plt.close()
    print "Created PDF plot file " + plotpath

# ----------------------------------------------------------------------
