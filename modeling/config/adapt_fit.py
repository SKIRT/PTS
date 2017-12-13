#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.modeling.decomposition.decomposition import scalelength_scaleheight_ratios, degeyter_ratio

# -----------------------------------------------------------------

# Load environment and model suite
environment = load_modeling_environment_cwd()
#suite = environment.static_model_suite

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# -----------------------------------------------------------------
