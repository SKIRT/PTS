#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.tools import filesystem as fs
from pts.modeling.fitting.component import get_run_names

# -----------------------------------------------------------------

# Set the modeling path
modeling_path = fs.cwd()

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# The fitting run for which to fit the SED
run_names = get_run_names(modeling_path)
if len(run_names) == 0: raise RuntimeError("No fitting runs found: first run configure_fit to create a new fitting run")
elif len(run_names) == 1: definition.add_fixed("name", "name of the fitting run", run_names[0])
else: definition.add_required("name", "string", "name of the fitting run", choices=run_names)

# Add optional arguments
definition.add_flag("visualise", "make visualisations")

# -----------------------------------------------------------------
