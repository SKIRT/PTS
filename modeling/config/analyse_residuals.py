#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.modeling.analysis.run import AnalysisRuns
from pts.modeling.core.environment import verify_modeling_cwd
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

modeling_path = verify_modeling_cwd()
runs = AnalysisRuns(modeling_path)

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# Positional optional
if runs.empty: raise ValueError("No analysis runs present (yet)")
elif runs.has_single: definition.add_fixed("run", "name of the analysis run", runs.single_name)
else: definition.add_positional_optional("run", "string", "name of the analysis run for which to launch the heating simulations", runs.last_name, runs.names)

# -----------------------------------------------------------------

definition.add_optional("nbins", "positive_integer", "number of bins for the distributions", 20)

# -----------------------------------------------------------------

# For clip mask
definition.add_optional("min_npixels", "positive_integer", "minimum number of pixels", 1)
definition.add_optional("connectivity", "positive_integer", "connectiviy", 4)

# -----------------------------------------------------------------

# REMAKE OR REPLOT
definition.add_flag("remake_residuals", "remake the residual maps")
definition.add_flag("remake_weighed", "remake the weighed residual maps")
definition.add_flag("remake_distributions", "remake the residuals distributions")
definition.add_flag("remake_weighed_distributions", "remake the weighed residuals distributions")
definition.add_flag("replot_distributions", "replot the residuals distributions")
definition.add_flag("replot_weighed_distributions", "replot the weighed residuals distributions")
definition.add_flag("replot_residuals", "replot residual maps")
definition.add_flag("replot_weighed", "replot weighed residual maps")

# -----------------------------------------------------------------

# For PNG
definition.add_optional("colours", "string", "colour or colour map for plotting", default=default_color)
definition.add_optional("scale", "string", "scaling", default_scale, scales)
definition.add_optional("interval", "string", "interval", default_interval)
definition.add_optional("alpha_method", "string", "alpha method", default_alpha_method, suggestions=alpha_methods)
definition.add_optional("peak_alpha", "real", "alpha of peak value", 1.5)

# -----------------------------------------------------------------
