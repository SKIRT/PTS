#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.eagle.observer Producing mock observations of the SKIRT results for a given SKIRT-run record
#
# The function in this module produces mock observations of the SKIRT results for a given SKIRT-runs database record.
# The mock observations include an SED plot, basic RGB images in the optical range, and (most importantly)
# an 'info' file containing fluxes for various bands in addition to other statistics.

# -----------------------------------------------------------------

import os
import os.path

from .skirtrun import SkirtRun
from ..core.plot.seds import plotseds
from ..core.plot.rgbimages import makergbimages
from .makeinfofile import makeinfofile

# -----------------------------------------------------------------

## This function produces mock observations of the SKIRT results for the given SKIRT-runs database record.
# The function uses the runid and snaptag fields in the specified database record.
# The mock observations include an SED plot, basic RGB images in the optical range, and (most importantly)
# an 'info' file containing fluxes for various bands in addition to other statistics.
#
# It is assumed that the appropriate SKIRT results directory contains the complete SKIRT simulation output.
# The files produced by this function are placed in the corresponding SKIRT-run visualization directory.
#
def observe(record):
    # get access to the appropriate SKIRT result directories
    skirtrun = SkirtRun(record["runid"])
    simulation = skirtrun.simulation()
    vispath = skirtrun.vispath()

    # create SED plot
    plotseds(simulation, output_path=vispath)

    # create basic RGB images at the SDSS gri wavelengths (not integrated over bands)
    makergbimages(simulation, wavelength_tuples=((0.753,0.617,0.470),), output_path=vispath)

    # get the redshift corresponding to the snapshot tag
    redshift = { 28:   0.0000,
                 27:   0.1006,
                 26:   0.1827,
                 25:   0.2709,
                 24:   0.3657,
                 23:   0.5031,
                 22:   0.6152,
                 21:   0.7356,
                 20:   0.8651,
                 19:   1.0041,
                 18:   1.2593,
                 17:   1.4867,
                 16:   1.7370,
                 15:   2.0124,
                 14:   2.2370,
                 13:   2.4784,
                 12:   3.0165,
                 11:   3.5280,
                 10:   3.9837,
                  9:   4.4852,
                  8:   5.0372,
                  7:   5.4874,
                  6:   5.9712,
                  5:   7.0496,
                  4:   8.0746,
                  3:   8.9879,
                  2:   9.9930,
                  1:  15.1323,
                  0:  20.0000 } [record["snaptag"]]

    # create the 'info' file
    makeinfofile(skirtrun, redshift=redshift)

# -----------------------------------------------------------------
