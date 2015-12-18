#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.makergbimages Make RGB images for the output of a SKIRT simulation.
#
# This script combines frames in each SKIRT \c prefix_instr_total.fits output file into
# an RGB image, placing a PNG file with the same name (but with the ".pdf" extension instead of ".fits")
# next to the original file. By default the script uses the first, last and middle frames in the fits files.
# To specify other frames/wavelengths you need to adjust the script.
#
# The script expects the complete output of a SKIRT simulation to be present (including log file etc.).
# If there are no arguments, the script processes all simulation output sets residing in the current directory.
# If the first argument contains a slash, the script processes all simulation output sets in the indicated directory.
# If the first argument does not contain a slash, the script processes just the simulation in the current directory
# with the indicated prefix.

# -----------------------------------------------------------------

# Import standard modules
import sys

# Import the relevant PTS classes and modules
from pts.core.simulation.simulation import createsimulations
from pts.core.plot.rgbimages import makergbimages

# -----------------------------------------------------------------

# the wavelengths for each of the RGB images to be created
wavelength_tuples = None                                    # use the default frames
#wavelength_tuples=((0.75,0.60,0.45),)                       # use approximate SDSS gri wavelengths
#wavelength_tuples = [ (0.77,0.55,0.33), (333,100,24) ]      # use these arbitrary wavelength sets

# -----------------------------------------------------------------

print "Starting makergbimages..."

# get the command-line argument specifying the simulation(s)
argument = sys.argv[1] if len(sys.argv) > 1 else ""

# construct the list of simulation objects and make the images
for simulation in createsimulations(argument):
    makergbimages(simulation, wavelength_tuples)

print "Finished makergbimages."

# -----------------------------------------------------------------
