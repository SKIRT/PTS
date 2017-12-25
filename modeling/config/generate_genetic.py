#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.core.environment import verify_modeling_cwd
from pts.modeling.config.component import definition

# -----------------------------------------------------------------

modeling_path = verify_modeling_cwd()

# -----------------------------------------------------------------

default_scale = "logarithmic"
scales = ["linear", "logarithmic"]

default_initial_generation_method = "random"
methods = ["random", "grid"]

# -----------------------------------------------------------------

definition = definition.copy()

# Add optional
definition.add_optional("ngenerations", "positive_integer", "number of generations to run in one run (ngenerations > 1 is only allowed for local execution)", 1)
definition.add_optional("nmodels", "positive_integer", "number of models per generation", 80)

definition.add_flag("manual_initial_generation", "generate intitial popultion manually", False)
definition.add_optional("manual_initial_generation_method", "string", "method for generating the initial generation manually", default=default_initial_generation_method, choices=methods)

# Scale
definition.add_optional("default_scale", "string", "default parameter scale (also for generating the initial generation manually)", default=default_scale, choices=scales)

# Check recurrence
definition.add_flag("check_recurrence", "check for recurrence of models that have been simulated previously", True)
definition.add_optional("recurrence_rtol", "positive_real", "relative tolerance for recurrence checking", 1e-5)
definition.add_optional("recurrence_atol", "positive_real", "absolute tolerance for recurrence checking", 1e-8)

# Make animations?
definition.add_flag("animate", "make animations", False)

# -----------------------------------------------------------------
