#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.eagle.plotgreybodyfit Plot a modified blackbody fit to
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

# import pts modules
from ..core.tools import archive as arch
from ..core.filter.broad import BroadBandFilter
from ..core.basics.greybody import GreyBody, kappa350_Cortese

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

    # load and plot the total SED
    filepath = simulation.seddatpaths()[0]
    lambdav, fluxv = np.loadtxt(arch.opentext(filepath), usecols=(0,1), unpack=True)
    lambdav = simulation.convert(lambdav, to_unit='micron', quantity='wavelength')
    fluxv = simulation.convert(fluxv, to_unit='Jy', quantity='fluxdensity', wavelength=lambdav)
    plot = (lambdav>=10) & (lambdav<=1000)
    plt.plot(lambdav[plot], fluxv[plot], color='b', label="SKIRT galaxy SED")

    # load and plot the contributions from HII particles (stellar emission) and gas particles (dust emission)
    # --> we do this inside a try block because these columns are not always available
    try:
        fstrdirv, fstrscav, ftotdusv = np.loadtxt(arch.opentext(filepath), usecols=(2,3,4), unpack=True)
        fstrdirv = simulation.convert(fstrdirv, to_unit='Jy', quantity='fluxdensity', wavelength=lambdav)
        fstrscav = simulation.convert(fstrscav, to_unit='Jy', quantity='fluxdensity', wavelength=lambdav)
        ftotdusv = simulation.convert(ftotdusv, to_unit='Jy', quantity='fluxdensity', wavelength=lambdav)
        plt.plot(lambdav[plot], fstrdirv[plot]+fstrscav[plot], color='c', ls="dashed", label="  contribution from HII regions")
        plt.plot(lambdav[plot], ftotdusv[plot], color='y', ls="dashed", label="  contribution from other dust")
    except:
        pass

    # load and plot the Herschel continuum data points (160, 250, 350, 500 micron)
    info = { }
    infofile = arch.listdir(skirtrun.vispath(), "_info.txt")[0]
    for line in arch.opentext(os.path.join(skirtrun.vispath(),infofile)):
        if not line.startswith("#"):
            key,dummy,value = line.split(None, 2)
            info[key] = float(value)
    waves = np.array( [ BroadBandFilter(fs).pivotwavelength() for fs in ("Pacs.red","SPIRE.PSW","SPIRE.PMW","SPIRE.PLW")] )
    fluxes = np.array(( info['instr_xy_fluxdensity_pacs_red_continuum'],
                        info['instr_xy_fluxdensity_spire_psw_continuum'],
                        info['instr_xy_fluxdensity_spire_pmw_continuum'],
                        info['instr_xy_fluxdensity_spire_plw_continuum']  ))
    sigmas = np.array(( 3,1,1,3 ))      # pacs is less sensitive; longer wavelength fluxes are harder to measure
    plt.scatter(waves, fluxes, color='r', marker='*', label="Mock PACS/SPIRE fluxes")

    # fit a grey body to the Herschel fluxes and plot the result
    greybody = GreyBody(simulation.instrumentdistance(), 2, kappa350_Cortese)
    T,M = greybody.fit(waves, fluxes, sigmas)
    plt.plot(lambdav[plot], greybody(lambdav[plot], T, M), color='m',
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
