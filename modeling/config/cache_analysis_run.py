#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.core.remote.host import all_host_ids
from pts.modeling.config.component import definition

# -----------------------------------------------------------------

# Load the modeling environment and analysis runs
environment = load_modeling_environment_cwd()
runs = environment.analysis_runs

# -----------------------------------------------------------------

definition = definition.copy()

# -----------------------------------------------------------------

# ANALYSIS RUN
if runs.empty: raise RuntimeError("No analysis runs are present (yet)")
elif runs.has_single: definition.add_fixed("run", "name of the analysis run", runs.single_name)
else: definition.add_positional_optional("run", "string", "name of the analysis run", runs.last_name, runs.names)

# Cache host ID
definition.add_positional_optional("remote", "string", "remote host ID to cache the analysis runs to", default=environment.cache_host_id, choices=all_host_ids())

# -----------------------------------------------------------------
