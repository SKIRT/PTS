#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.plotprogress Plot progress for the various phases of a SKIRT simulation.
#
# This script plots the progress in function of time for certain phases of a SKIRT simulation, based on the log
# messages. A seperate PDF plot is created for each of the following phases, if present in the simulation:
# - shooting photons for stellar emission ("prefix_progress_stellar_photons.pdf");
# - calculating dust emission spectra ("prefix_progress_dust_spectra.pdf");
# - shooting photons for dust emission ("prefix_progress_dust_photons.pdf");
#
# The dust self-absorption phase, if present, is ignored in the current implementation of the script.
#
# For multi-process (MPI) simulations with verbose logging (i.e. with a separate log file per process),
# the progress for all processes is displayed on the same plot.
#
# The script expects the complete output of a SKIRT simulation to be present (including log file etc.).
# If there are no arguments, the script processes all simulation output sets residing in the current directory.
# If the first argument contains a slash, the script processes all simulation output sets in the indicated directory.
# If the first argument does not contain a slash, the script processes just the simulation in the current directory
# with the indicated prefix.
#

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os

# Import the relevant PTS classes and modules
from pts.core.simulation.simulation import createsimulations
from pts.core.extract.progress import ProgressExtractor
from pts.core.plot.progress import ProgressPlotter

# -----------------------------------------------------------------

# Look for a file in the current working directory that contains extracted progress information
progress_file_path = os.path.join(os.getcwd(), "progress.dat")
if os.path.isfile(progress_file_path):

    # Create a ProgressExtractor instance from the saved progress data
    extractor = ProgressExtractor.open_table(progress_file_path)

# If extracted timeline information is not present, first perform the extraction
else:

    # Create a SkirtSimulation object based on a log file present in the current working directory
    simulation = createsimulations(single=True)

    # Create a new ProgressExtractor instance
    extractor = ProgressExtractor()

    # Run the extractor
    extractor.run(simulation)

# Determine the path to the plot file
plot_path = os.path.join(os.getcwd(), "progress.pdf")

# Create a ProgressPlotter instance
plotter = ProgressPlotter()

# Run the progress plotter
plotter.run(extractor.table, plot_path)

# -----------------------------------------------------------------
