#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.core.environment import verify_modeling_cwd
from pts.modeling.analysis.run import AnalysisRuns
from pts.magic.core.rgba import alpha_methods
from pts.modeling.config.component import definition

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

definition = definition.copy()

# Positional optional
if runs.empty: raise ValueError("No analysis runs present (yet)")
elif runs.has_single: definition.add_fixed("run", "name of the analysis run", runs.single_name)
else: definition.add_positional_optional("run", "string", "name of the analysis run", runs.last_name, runs.names)

# -----------------------------------------------------------------

# For reference SED
definition.add_flag("not_clipped", "use the observed fluxes from the truncated (not clipped) images")
definition.add_optional("additional_error", "percentage", "additional percentual error for the observed flux points", 0.1)

# -----------------------------------------------------------------

definition.add_optional("nbins", "positive_integer", "number of bins for the distributions", 20)

# -----------------------------------------------------------------

# For clip mask
definition.add_optional("min_npixels", "positive_integer", "minimum number of pixels", 1)
definition.add_optional("connectivity", "positive_integer", "connectiviy", 4)

# -----------------------------------------------------------------

# For PNG
definition.add_optional("colours", "string", "colour or colour map for plotting", default=default_color)
definition.add_optional("scale", "string", "scaling", default_scale, scales)
definition.add_optional("interval", "string", "interval", default_interval)
definition.add_optional("alpha_method", "string", "alpha method", default_alpha_method, suggestions=alpha_methods)
definition.add_optional("peak_alpha", "real", "alpha of peak value", 1.5)

# -----------------------------------------------------------------

# Plot?
definition.add_flag("plot", "do plotting", True)

# -----------------------------------------------------------------
