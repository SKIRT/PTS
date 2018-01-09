#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.modeling.config.component import definition

# -----------------------------------------------------------------

# Load environment and model suite
environment = load_modeling_environment_cwd()
runs = environment.fitting_runs

# -----------------------------------------------------------------

definition = definition.copy()

# -----------------------------------------------------------------

# Name for the reweighing
definition.add_required("name", "string", "name for the reweiging")

# The fitting run for which to adapt the configuration
if runs.empty: raise RuntimeError("No fitting runs are present")
elif runs.has_single: definition.add_fixed("fitting_run", "name of the fitting run", runs.single_name)
else: definition.add_required("fitting_run", "string", "name of the fitting run", choices=runs.names)

# Generations?

# -----------------------------------------------------------------

definition.add_optional("filters", "filter_list", "filters to use for the evaluation (None means default fitting filters)")

# -----------------------------------------------------------------
