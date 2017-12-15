#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.remove_generations Remove certain generations from a fitting run.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import verify_modeling_cwd
from pts.modeling.fitting.run import FittingRuns
from pts.core.tools import filesystem as fs
from pts.core.basics.log import log

# -----------------------------------------------------------------

# Get the modeling commands
modeling_path = verify_modeling_cwd()
runs = FittingRuns(modeling_path)

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# The fitting run for which to explore the parameter space
if runs.empty: raise RuntimeError("No fitting runs are present")
elif runs.has_single: definition.add_fixed("name", "name of the fitting run", runs.single_name)
else: definition.add_required("name", "string", "name of the fitting run", choices=runs.names)

# Generations to remove
definition.add_positional_optional("generations", "string_list", "generations to remove (none means all)")

# Get configuration
config = parse_arguments("remove_generations", definition)

# -----------------------------------------------------------------

# Load the fitting run
fitting_run = runs.load(config.name)

# Loop over the generations
generation_names = config.generations if config.generations is not None else fitting_run.generation_names
for generation_name in generation_names:

    # Inform the user
    log.info("Removing generation '" + generation_name + "' ...")

    # Remove this generation from the generations table
    fitting_run.generations_table.remove_entry(generation_name)
    fitting_run.generations_table.save()

    # Get generation path
    generation_path = fitting_run.get_generation_path(generation_name)

    # Remove the generation directory
    fs.remove_directory(generation_path)

# -----------------------------------------------------------------
