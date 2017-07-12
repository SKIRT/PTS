#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.modeling.fitting.run import has_single_fitting_run, get_fitting_run_names, has_fitting_runs, get_single_fitting_run_name
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

# Set the modeling path
modeling_path = fs.cwd()

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# -----------------------------------------------------------------

# No fitting runs (yet)
if not has_fitting_runs(modeling_path):

    # Check whether we can prompt for model definitions


# One fitting run
elif has_single_fitting_run(modeling_path): definition.add_fixed("fitting_run", "string", get_single_fitting_run_name)

# Multiple fitting runs
else: definition.add_required("fitting_run", "string", "name of the fitting run to use", choices=get_fitting_run_names(modeling_path))

# -----------------------------------------------------------------

# Add flags
definition.add_flag("show", "show", False)
definition.add_flag("replot", "make plots again", True)

# -----------------------------------------------------------------
