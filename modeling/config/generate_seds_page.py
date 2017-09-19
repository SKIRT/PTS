#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.config.generate_page import definition
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.magic.core.rgba import alpha_methods

# -----------------------------------------------------------------

default_alpha_method = "combined"
default_scale = "log"
default_color = "jet"
default_mask_color = "black"

scales = ["log", "sqrt"]
default_colour = "jet"
default_interval = "pts"

# -----------------------------------------------------------------

environment = load_modeling_environment_cwd()
analysis_runs = environment.analysis_runs

# -----------------------------------------------------------------

definition.add_required("analysis_run", "string", "analysis run to clone", choices=analysis_runs.names)

# -----------------------------------------------------------------

# For PNG
definition.add_optional("colours", "string", "colour or colour map for plotting", default=default_color)
definition.add_optional("scale", "string", "scaling", default_scale, scales)
definition.add_optional("interval", "string", "interval", default_interval)
definition.add_optional("alpha_method", "string", "alpha method", default_alpha_method, suggestions=alpha_methods)
definition.add_optional("peak_alpha", "real", "alpha of peak value", 1.5)

# -----------------------------------------------------------------
