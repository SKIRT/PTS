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

    # create the 'info' file
    makeinfofile(skirtrun, record["snaptag"])

# -----------------------------------------------------------------
