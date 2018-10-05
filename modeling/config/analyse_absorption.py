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

definition.add_flag("show", "showing", True)
definition.add_flag("plot", "plotting", True)

# -----------------------------------------------------------------

# Recalculate properties
definition.add_flag("recalculate", "recalculate the properties")
definition.add_flag("recalculate_total", "recalculate total properties")
definition.add_flag("recalculate_bulge", "recalculate bulge properties")
definition.add_flag("recalculate_disk", "recalculate disk properties")
definition.add_flag("recalculate_old", "recalculate old properties")
definition.add_flag("recalculate_young", "recalculate young properties")
definition.add_flag("recalculate_sfr", "recalculate SFR properties")
definition.add_flag("recalculate_unevolved", "recalculate unevolved properties")

# Replot
definition.add_flag("replot", "remake the plots")
definition.add_flag("replot_total", "replot total")
definition.add_flag("replot_bulge", "replot bulge")
definition.add_flag("replot_disk", "replot disk")
definition.add_flag("replot_old", "replot old")
definition.add_flag("replot_young", "replot young")
definition.add_flag("replot_sfr", "replot SFR")
definition.add_flag("replot_unevolved", "replot unevolved")

# -----------------------------------------------------------------

definition.add_optional("unit", "photometric_unit", "SED plotting unit", "W/m2 [density]", convert_default=True)

# -----------------------------------------------------------------
