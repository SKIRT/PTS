#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.fitting.run import has_single_fitting_run, get_fitting_run_names, has_fitting_runs, get_single_fitting_run_name
from pts.core.tools import filesystem as fs
from pts.modeling.config.generate_page import definition

# -----------------------------------------------------------------

# Set the modeling path
modeling_path = fs.cwd()

# -----------------------------------------------------------------

# Fitting run setting
if not has_fitting_runs(modeling_path): raise RuntimeError("There are no fitting runs")
elif has_single_fitting_run(modeling_path): definition.add_fixed("fitting_run", "string", get_single_fitting_run_name)
else: definition.add_required("fitting_run", "string", "name of the fitting run to use", choices=get_fitting_run_names(modeling_path))

# -----------------------------------------------------------------

# Flags
definition.add_flag("show", "show the page", False)

# -----------------------------------------------------------------
