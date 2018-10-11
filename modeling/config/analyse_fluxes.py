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

# Positional optional
if runs.empty: raise ValueError("No analysis runs present (yet)")
elif runs.has_single: definition.add_fixed("run", "name of the analysis run", runs.single_name)
else: definition.add_positional_optional("run", "string", "name of the analysis run", runs.last_name, runs.names)

# -----------------------------------------------------------------

# Unit for mock fluxes
definition.add_optional("unit", "photometric_unit", "unit to use for the fluxes", "Jy", convert_default=True)

# -----------------------------------------------------------------

# Plot?
definition.add_flag("plot", "do plotting", True)

# -----------------------------------------------------------------

# Additional relative error for observed SED
definition.add_optional("additional_error", "percentage", "additional percentual error for the observed flux points", "10", convert_default=True) # just having 0.1 as default DOES NOT WORK (will be converted into 0.001): THIS IS AN OPEN BUG IN PTS

# -----------------------------------------------------------------
