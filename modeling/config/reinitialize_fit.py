#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.core.environment import verify_modeling_cwd
from pts.modeling.fitting.run import FittingRuns
from pts.modeling.build.suite import ModelSuite
from pts.modeling.component.component import get_default_fitting_method
from pts.modeling.config.component import definition

# -----------------------------------------------------------------

# Determine the modeling path
modeling_path = verify_modeling_cwd()
suite = ModelSuite.from_modeling_path(modeling_path)
runs = FittingRuns(modeling_path)

# -----------------------------------------------------------------

origins = ["model", "fitting_run", "analysis_run"]

# Fitting methods
default_fitting_method = get_default_fitting_method(modeling_path) # as defined in the modeling configuration
fitting_methods = ["genetic", "grid"]

# -----------------------------------------------------------------

definition = definition.copy()

# FITTING RUN
if runs.empty: definition.add_required("name", "string", "name for the fitting run")
else: definition.add_required("name", "string", "name for the fitting run", forbidden=runs.names)

# -----------------------------------------------------------------

definition.add_required("origin", "string", "origin of the fitting run", choices=origins)

# -----------------------------------------------------------------

# FITTING METHOD
definition.add_optional("fitting_method", "string", "fitting method", default_fitting_method, fitting_methods)

# -----------------------------------------------------------------

# Sections
#definition.add_section("ranges", "parameter ranges")

# -----------------------------------------------------------------
