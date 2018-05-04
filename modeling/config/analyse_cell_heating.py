#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.core.environment import verify_modeling_cwd
from pts.modeling.analysis.run import AnalysisRuns
from pts.modeling.config.component import definition

# -----------------------------------------------------------------

modeling_path = verify_modeling_cwd()
runs = AnalysisRuns(modeling_path)

# -----------------------------------------------------------------

definition = definition.copy()

# ANALYSIS RUN
if runs.empty: raise ValueError("No analysis runs present (yet)")
elif runs.has_single: definition.add_fixed("run", "name of the analysis run", runs.single_name)
else: definition.add_positional_optional("run", "string", "name of the analysis run for which to launch the heating simulations", runs.names)

# -----------------------------------------------------------------

# The number of bins
definition.add_optional("nbins", "positive_integer", "number of bins", 20)
definition.add_optional("nradial_bins", "positive_integer", "number of radial bins", 200)

# -----------------------------------------------------------------

# Plot
definition.add_flag("plot", "do plotting", True)
definition.add_flag("plot_distribution", "plot distribution", True)
definition.add_flag("plot_radial_distribution", "plot radial distribution", True)
definition.add_flag("plot_map", "plot map", True)

# -----------------------------------------------------------------

# Recreate
definition.add_flag("recreate_table", "recreate the absorption table")
definition.add_flag("recalculate_fractions", "recalculate the heating fractions")
definition.add_flag("recalculate_distributions", "recreate the heating distributions")
definition.add_flag("recalculate_distribution", "recalculate")
definition.add_flag("recalculate_distribution_diffuse", "recalculate")
definition.add_flag("recalculate_radial_distribution", "recalculate")
definition.add_flag("replot", "replot")
definition.add_flag("replot_distribution", "replot the distribution")
definition.add_flag("replot_radial_distribution", "replot the radial distribution")
definition.add_flag("replot_map", "replot the map")

# -----------------------------------------------------------------
