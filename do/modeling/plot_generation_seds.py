#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.plot_generation_seds Plot the seds of all simulations of a generation.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.modeling.fitting.statistics import FittingStatistics

# -----------------------------------------------------------------

environment = load_modeling_environment_cwd()
runs = environment.fitting_runs

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# The fitting run name
if runs.empty: raise ValueError("No fitting runs present (yet)")
elif runs.has_single: definition.add_fixed("run", "name of the fitting run", runs.single_name)
else: definition.add_required("run", "string", "name of the fitting run", choices=runs.names)

# Generation
definition.add_required("generation", "string", "generation name")

# Number of random simulations
definition.add_optional("random", "positive_integer", "pick a specified number of random simulations to plot")

# Additional relative error
definition.add_optional("additional_error", "percentage", "additional percentual error for the observed flux points")

# Output file
definition.add_optional("output", "string", "output file name")

# Get configuration
config = parse_arguments("plot_generation_seds", definition, "Plot the seds of all simulations of a generation")

# -----------------------------------------------------------------

# Set plot command
plot_command = 'plot seds "' + config.generation + '"'
if config.random is not None: plot_command += " --random " + str(config.random)
if config.additional_error is not None: plot_command += " --additional_error " + str(config.additional_error * 100)
if config.output is not None: plot_command += ' --path "' + config.output + '"'

# -----------------------------------------------------------------

# Create the fitting statistics
statistics = FittingStatistics()

# Set properties
statistics.config.run = config.run
statistics.config.interactive = False
statistics.config.commands = [plot_command]

# Run
statistics.run()

# -----------------------------------------------------------------
