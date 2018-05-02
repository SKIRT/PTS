#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.core.environment import verify_modeling_cwd
from pts.modeling.fitting.run import FittingRuns
from pts.modeling.config.component import definition

# -----------------------------------------------------------------

# Set the modeling path
modeling_path = verify_modeling_cwd()
runs = FittingRuns(modeling_path)

# -----------------------------------------------------------------

definition = definition.copy()

# The fitting run for which to fit the SED
if runs.empty: raise RuntimeError("No fitting runs are present (yet)")
elif runs.has_single: definition.add_fixed("name", "name of the fitting run", runs.single_name)
else: definition.add_required("name", "string", "name of the fitting run", choices=runs.names)

# Add optional arguments
definition.add_flag("visualise", "make visualisations")

# -----------------------------------------------------------------

# Flags
definition.add_flag("unfinished", "also include unfinished generations")
definition.add_optional("rerun", "string", "rerun for a certain generation")
definition.add_flag("rerun_all", "rerun all")

# -----------------------------------------------------------------

# Flags
definition.add_flag("write", "write", True)
definition.add_flag("plot", "make plots", True)

# -----------------------------------------------------------------

definition.add_flag("per_generation", "only calculate probabilities per generation (because e.g. the fitting circumstances have changed)")

# -----------------------------------------------------------------

# Re-calculate the chi squared values
definition.add_flag("recalculate_chisquared", "recalculate the chi squared values", False)
definition.add_optional("recalculate_generations", "string_list", "recalculate the chi squared values only for these generations")
definition.add_optional("recalculate_simulations", "string_list", "recalculate the chi squared values only for these simulations")

# -----------------------------------------------------------------
