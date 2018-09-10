#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import standard modules
import warnings

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
if runs.empty: warnings.warn("No analysis runs present (yet)")
elif runs.has_single: definition.add_fixed("run", "name of the analysis run", runs.single_name)
else: definition.add_positional_optional("run", "string", "name of the analysis run", runs.names)

# -----------------------------------------------------------------

# Plotting
definition.add_flag("plot", "do plotting", True)
definition.add_flag("replot", "replot all")
definition.add_flag("replot_sfr", "replot all SFR maps")
definition.add_flag("replot_mass", "replot all mass maps")
definition.add_flag("replot_ssfr", "replot all sSFR maps")

# Projected
definition.add_flag("replot_projected", "replot projected", False)
definition.add_flag("replot_projected_sfr", "replot projected SFR", False)
definition.add_flag("replot_projected_mass", "replot projected mass", False)
definition.add_flag("replot_projected_ssfr", "replot projected sSFR", False)

# Cell maps
definition.add_flag("replot_cell_maps", "replot cell maps", False)
definition.add_flag("replot_cell_ssfr_maps", "replot cell sSFR maps", False)

# -----------------------------------------------------------------

definition.add_flag("project", "create and plot maps from the 3D data", True)

# -----------------------------------------------------------------
