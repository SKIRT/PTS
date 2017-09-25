#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.show_analysis_runs Show the analysis runs.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import load_modeling_environment_cwd

# -----------------------------------------------------------------

# Determine the modeling path
environment = load_modeling_environment_cwd()
#runs = environment.analysis_runs
context = environment.analysis_context

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Get configuration
config = parse_arguments("show_model", definition)

# -----------------------------------------------------------------

# Loop over the cache host IDS
for host_id in context.cache_host_ids:

    print("")
    print(host_id.upper() + ":")
    print("")

    # Get the run names
    run_names = context.get_run_names_for_host_id(host_id)

    for name in run_names:
        print(" - " + name)

# -----------------------------------------------------------------

# Show the local analysis runs

print("")
print("LOCAL:")
print("")

for name in context.analysis_run_names:
    print(" - " + name)

print("")

# -----------------------------------------------------------------
