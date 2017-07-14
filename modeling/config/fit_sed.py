#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.modeling.core.environment import verify_modeling_cwd
from pts.modeling.fitting.run import FittingRuns

# -----------------------------------------------------------------

# Set the modeling path
modeling_path = verify_modeling_cwd()
runs = FittingRuns(modeling_path)

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# The fitting run for which to fit the SED
# FITTING RUN
if runs.empty: raise RuntimeError("No fitting runs are present (yet)")
elif runs.has_single: definition.add_fixed("name", "name of the fitting run", runs.single_name)
else: definition.add_required("name", "string", "name of the fitting run", choices=runs.names)

# Add optional arguments
definition.add_flag("visualise", "make visualisations")

# -----------------------------------------------------------------
