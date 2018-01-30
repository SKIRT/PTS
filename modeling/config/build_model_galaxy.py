#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.modeling.config.component import definition

# -----------------------------------------------------------------

# Load environment and model suite
environment = load_modeling_environment_cwd()
suite = environment.static_model_suite

# -----------------------------------------------------------------

definition = definition.copy()

# Add settings
if suite.has_models: definition.add_required("name", "string", "name for the model", forbidden=suite.model_names)
else: definition.add_positional_optional("name", "string", "name for the model", default="standard")

# Description
definition.add_positional_optional("description", "string", "description of the model")

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

definition.add_flag("use_defaults", "use default parameter values", False)

# -----------------------------------------------------------------

definition.add_flag("prompt_parameters", "fill in the global SED fit parameters manually (don't try to fetch from DustPedia)", False)

# -----------------------------------------------------------------
