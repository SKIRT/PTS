#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.modeling.core.environment import verify_modeling_cwd

# -----------------------------------------------------------------

modeling_path = verify_modeling_cwd()

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# -----------------------------------------------------------------

#origins = ["initial_model", "fitting_run", "custom_model"]
origins = ["model", "fitting_run"]

# -----------------------------------------------------------------

definition.add_required("origin", "string", "origin of the analysis model", choices=origins)

# -----------------------------------------------------------------
