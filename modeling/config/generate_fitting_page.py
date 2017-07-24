#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.config.generate_page import definition
from pts.modeling.fitting.run import FittingRuns
from pts.modeling.core.environment import verify_modeling_cwd

# -----------------------------------------------------------------

# Set the modeling path
modeling_path = verify_modeling_cwd()
runs = FittingRuns(modeling_path)

# -----------------------------------------------------------------

# Fitting run setting
if runs.empty: raise RuntimeError("No fitting runs are present (yet)")
elif runs.has_single: definition.add_fixed("fitting_run", "name of the fitting run", runs.single_name)
else: definition.add_required("fitting_run", "string", "name of the fitting run", choices=runs.names)

# -----------------------------------------------------------------
