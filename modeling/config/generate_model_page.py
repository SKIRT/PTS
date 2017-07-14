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
if not has_fitting_runs(modeling_path): raise RuntimeError("There are no fitting runs")
elif has_single_fitting_run(modeling_path): definition.add_fixed("fitting_run", "string", get_single_fitting_run_name)
else: definition.add_required("fitting_run", "string", "name of the fitting run to use", choices=get_fitting_run_names(modeling_path))

# -----------------------------------------------------------------

# Flags
definition.add_flag("show", "show the page", False)

# -----------------------------------------------------------------
