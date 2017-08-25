#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.evolve.show_reproductions Show the reproductions.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import GalaxyModelingEnvironment
from pts.modeling.fitting.component import FittingRun
from pts.core.tools import formatting as fmt
from pts.modeling.core.environment import verify_modeling_cwd
from pts.modeling.fitting.run import FittingRuns

# -----------------------------------------------------------------

modeling_path = verify_modeling_cwd()
runs = FittingRuns(modeling_path)

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()

# -----------------------------------------------------------------

# THE FITTING RUN
if runs.empty: raise RuntimeError("No fitting runs are present (yet)")
elif runs.has_single: definition.add_fixed("fitting_run", "name of the fitting run", runs.single_name)
else: definition.add_required("fitting_run", "string", "name of the fitting run", choices=runs.names)

# -----------------------------------------------------------------

# Generation
definition.add_positional_optional("generations", "string_list", "name of the generations for which to show the reproductions")

# Create the configuration
config = parse_arguments("show_reproductions", definition)

# -----------------------------------------------------------------

# Load the modeling environment
environment = GalaxyModelingEnvironment(modeling_path)

# -----------------------------------------------------------------

# Load the fitting run
fitting_run = FittingRun.from_name(modeling_path, config.fitting_run)
#print(fitting_run.parameter_base_types)

# -----------------------------------------------------------------

# Get generation names
generations = config.generations if config.generations is not None else fitting_run.genetic_generations

# -----------------------------------------------------------------

print("")

# Loop over the generations
for generation_name in generations:

    # Print generation name
    print(fmt.bold + fmt.cyan + generation_name + fmt.reset)
    print("")

    # Get the generation platform
    platform = fitting_run.get_generation_platform(generation_name)

    # -----------------------------------------------------------------

    # Show
    platform.show_reproductions()

# -----------------------------------------------------------------
