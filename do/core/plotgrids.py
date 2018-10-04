#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.plotgrids Plot dust grids for SKIRT output files.
#
# This script plots dust grids for SKIRT \c prefix_ds_gridxx.dat output files, placing PDF files with the same
# name (but with the ".pdf" extension instead of ".dat") next to the original files.
#
# The script expects the complete output of a SKIRT simulation to be present (including log file etc.).
# If there are no arguments, the script processes all simulation output sets residing in the current directory.
# If the first argument contains a slash, the script processes all simulation output sets in the indicated directory.
# If the first argument does not contain a slash, the script processes just the simulation in the current directory
# with the indicated prefix.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.simulation.simulation import createsimulations
from pts.core.plot.grids import plotgrids
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()
definition.add_positional_optional("simulation", "string", "simulation specification", default="")
definition.add_optional("linewidth", "positive_real", "line width", 0.1)
definition.add_optional("maxlevel", "positive_integer", "maximum tree level")

# Read the command line arguments
config = parse_arguments("plotgrids", definition, description="Unmount a remote mounted with PTS")

# -----------------------------------------------------------------

print("Starting plotgrids...")

# construct the list of simulation objects and make the plots
for simulation in createsimulations(config.simulation): plotgrids(simulation, linewidth=config.linewidth, maxlevel=config.maxlevel)

print("Finished plotgrids")

# -----------------------------------------------------------------
