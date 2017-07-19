#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.config.generate_page import definition
from pts.modeling.analysis.run import AnalysisRuns
from pts.modeling.core.environment import verify_modeling_cwd

# -----------------------------------------------------------------

# Set the modeling path
modeling_path = verify_modeling_cwd()
runs = AnalysisRuns(modeling_path)

# -----------------------------------------------------------------

# The analysis run
if runs.empty: raise ValueError("There are no analysis runs (yet)")
elif runs.has_single: definition.add_fixed("analysis_run", "analysis run", runs.single_name)
else: definition.add_required("analysis_run", "string", "analysis run", choices=runs.names)

# -----------------------------------------------------------------
