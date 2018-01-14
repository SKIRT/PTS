#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.plot_fitting_filters Plot the fitting filters of a fitting run.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.core.plot.transmission import TransmissionPlotter

# -----------------------------------------------------------------

environment = load_modeling_environment_cwd()
runs = environment.fitting_runs

# -----------------------------------------------------------------

# Create definition
definition = ConfigurationDefinition()

# FITTING RUN
if runs.empty: raise RuntimeError("No fitting runs are present (yet)")
elif runs.has_single: definition.add_fixed("fitting_run", "name of the fitting run", runs.single_name)
else: definition.add_required("fitting_run", "string", "name of the fitting run", choices=runs.names)

# Get the arguments
config = parse_arguments("plot_fitting_filters", definition, "Plot the fitting filters of a fitting run")

# -----------------------------------------------------------------

# Load the fitting run
fitting_run = runs.load(config.fitting_run)

# -----------------------------------------------------------------

# Initialize the plotter
plotter = TransmissionPlotter()

# Add the filters
plotter.add_filters(fitting_run.fitting_filters)

# Run the plotter
plotter.run()

# -----------------------------------------------------------------
