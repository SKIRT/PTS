#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.evolve.generations Show generations.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.tools import formatting as fmt
from pts.modeling.core.environment import verify_modeling_cwd
from pts.modeling.fitting.run import FittingRuns

# -----------------------------------------------------------------

modeling_path = verify_modeling_cwd()
runs = FittingRuns(modeling_path)

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()

# FITTING RUN
if runs.empty: raise RuntimeError("No fitting runs are present (yet)")
else: definition.add_positional_optional("runs", "string_list", "names of the fitting runs", runs.names, choices=runs.names)

# Create the configuration
config = parse_arguments("generations", definition)

# -----------------------------------------------------------------

# Loop over the fitting runs
for name in config.runs:

    # Get the fitting run
    fitting_run = runs.load(name)

    # Show
    print("")
    print(fmt.green + fmt.underlined + name + fmt.reset + ":")
    print("")

    # Loop over the generations
    for generation_name in fitting_run.generation_names:
        print(" - " + generation_name)

    print("")

# -----------------------------------------------------------------
