#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.plotprogress Plot progress for the various phases of a SKIRT simulation.
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

# Import the relevant PTS classes and modules
from pts.core.extract.progress import extract_progress_cwd, ProgressTable
from pts.core.plot.progress import ProgressPlotter
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add flags
definition.add_flag("table", "save the extracted progress table")

# Get configuration
config = parse_arguments("plotprogress", definition)

# -----------------------------------------------------------------

# Look for a file in the current working directory that contains extracted progress information
progress_table_path = fs.join(fs.cwd(), "progress.dat")
if fs.is_file(progress_table_path): table = ProgressTable.from_file(progress_table_path)

# If extracted progress information is not present, first perform the extraction
else: table = extract_progress_cwd()

# -----------------------------------------------------------------

if config.table and not fs.is_file(progress_table_path): table.saveto(progress_table_path)

# -----------------------------------------------------------------

# Determine the path to the plotting directory
plot_path = fs.join(fs.cwd())

# Create a ProgressPlotter instance
plotter = ProgressPlotter()

# Run the progress plotter
plotter.run(table, plot_path)

# -----------------------------------------------------------------
