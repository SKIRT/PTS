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

definition.add_flag("topcat", "make plots using topcat (STILTS)")

# -----------------------------------------------------------------

# Replot
definition.add_flag("replot", "remake the plots")
definition.add_flag("replot_ssfr_funev", "remake the sSFR-Funev correlation plots")
definition.add_flag("replot_ssfr_salim_funev", "Salim")
definition.add_flag("replot_ssfr_ke_funev", "K&E")
definition.add_flag("replot_ssfr_mappings_ke_funev", "MAPPINGS + K&E")

# -----------------------------------------------------------------
