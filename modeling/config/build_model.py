#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.modeling.build.suite import ModelSuite
from pts.core.tools import filesystem as fs
from pts.modeling.config.build_stars import default_sfr
from pts.modeling.config.build_dust import default_dust_mass

# -----------------------------------------------------------------

modeling_path = fs.cwd()
suite = ModelSuite.from_modeling_path(modeling_path)

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# Add settings
model_names = suite.model_names
if len(model_names) == 0: definition.add_positional_optional("name", "string", "name for the model", default="standard")
else: definition.add_required("name", "string", "name for the model", forbidden=model_names)

# Description
definition.add_positional_optional("description", "string", "description of the model")

# Settings for stellar and dust components
definition.add_optional("sfr", "real", "average star formation rate", default_sfr)
definition.add_optional("dust_mass", "quantity", "estimated mass of the dust disk", default_dust_mass)

# -----------------------------------------------------------------
