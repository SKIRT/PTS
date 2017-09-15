#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.modeling.config.build_stars import default_sfr
from pts.modeling.config.build_dust import default_dust_mass
from pts.modeling.core.environment import load_modeling_environment_cwd

# -----------------------------------------------------------------

# Load environment and model suite
environment = load_modeling_environment_cwd()
suite = environment.static_model_suite

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# Add settings
if suite.has_models: definition.add_required("name", "string", "name for the model", forbidden=suite.model_names)
else: definition.add_positional_optional("name", "string", "name for the model", default="standard")

# Description
definition.add_positional_optional("description", "string", "description of the model")

# Settings for stellar and dust components
definition.add_optional("sfr", "real", "average star formation rate", default_sfr)
definition.add_optional("dust_mass", "quantity", "estimated mass of the dust disk", default_dust_mass)

# Ask for additional components
definition.add_flag("additional", "ask for additional components", True)

# -----------------------------------------------------------------

# ADVANCED: create model from a previous model
if suite.has_models: definition.add_optional("from_previous", "string", "create from previous model", choices=suite.model_names)
else: definition.add_fixed("from_previous", "create from previous model", None)

# -----------------------------------------------------------------

# Force overwrite
definition.add_flag("overwrite", "overwrite possibly existing model with this name", False)

# Show after adjusting/building
definition.add_flag("show", "show the components after the model is built", True)

# -----------------------------------------------------------------
