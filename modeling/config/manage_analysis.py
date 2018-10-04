#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.core.config.manage_simulations import definition

# -----------------------------------------------------------------

# Get the analysis runs
environment = load_modeling_environment_cwd()
runs = environment.analysis_runs

# -----------------------------------------------------------------

definition = definition.copy()

# -----------------------------------------------------------------

# THE ANALYSIS RUN
if runs.empty: raise ValueError("No analysis runs present (yet)")
elif runs.has_single: definition.add_fixed("run", "name of the analysis run", runs.single_name)
else: definition.add_required("run", "string", "name of the analysis run", runs.names)

# -----------------------------------------------------------------
