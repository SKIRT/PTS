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

# -----------------------------------------------------------------

# The analysis run
if runs.empty: print("No analysis runs present (yet)")
elif runs.has_single: definition.add_fixed("run", "name of the analysis run", runs.single_name)
else: definition.add_positional_optional("run", "string", "name of the analysis run", runs.last_name, runs.names)

# -----------------------------------------------------------------

#definition.add_optional("scale_heights", "real", "number of times to take the old stellar scale height as the vertical radius of the model", 15.)

# -----------------------------------------------------------------

definition.add_flag("earth", "write earth maps", True)
definition.add_flag("faceon", "write faceon maps", True)
definition.add_flag("edgeon", "write edgeon maps", True)

# -----------------------------------------------------------------

# Flags
definition.add_flag("plot", "do plotting", True)

# -----------------------------------------------------------------
