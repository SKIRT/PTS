#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.plot_with_reference_seds Plot a certain simulated SED with the modeling reference SED.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.modeling.core.environment import find_modeling_environment_up_cwd
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.plot.sed import SEDPlotter
from pts.core.data.sed import SED, ObservedSED, is_from_skirt
from pts.core.basics.plot import mpl, plotting_libraries

# -----------------------------------------------------------------

environment = find_modeling_environment_up_cwd()

# -----------------------------------------------------------------

# Create the definition
definition = ConfigurationDefinition()
definition.add_required("sed", "file_path", "path to the sed file")
definition.add_positional_optional("outfile_path", "string", "output file path")
definition.add_optional("library", "string", "plotting library", mpl, choices=plotting_libraries)

# Get the confguration
config = parse_arguments("plot_with_reference_sed", definition)

# -----------------------------------------------------------------

# Initialize the plotter
plotter = SEDPlotter()
plotter.config.library = config.library

# -----------------------------------------------------------------

# Load the modeled SED
if is_from_skirt(config.sed):

    sed = SED.from_skirt(config.sed)
    label = "Simulation"

# Load the mock observed SED
else:

    sed = ObservedSED.from_file(config.sed)
    label = "Mock observation"

# -----------------------------------------------------------------

# Add the SEDS
plotter.add_sed(environment.observed_sed, "Observation")
plotter.add_sed(sed, label)

# -----------------------------------------------------------------

# Plot
plotter.run(output=config.outfile_path)

# -----------------------------------------------------------------
