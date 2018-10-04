#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.analysis.run import AnalysisRuns
from pts.modeling.core.environment import verify_modeling_cwd
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

# Remake?
definition.add_flag("remake", "remake already existing maps", False)

# -----------------------------------------------------------------

definition.add_optional("hot_factor_range", "real_range", "range of factor to create the hot dust maps", "0.2>0.7", convert_default=True)
definition.add_optional("factor_nvalues", "positive_integer", "number of factors", 8)

# -----------------------------------------------------------------
