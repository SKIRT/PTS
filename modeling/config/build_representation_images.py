#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.modeling.build.suite import ModelSuite
from pts.modeling.core.environment import verify_modeling_cwd
from pts.modeling.config.component import definition

# -----------------------------------------------------------------

# Determine the modeling path
modeling_path = verify_modeling_cwd()
suite = ModelSuite.from_modeling_path(modeling_path)

# -----------------------------------------------------------------

definition = definition.copy()

# Name of the representation
representation_names = suite.representation_names
if len(representation_names) == 0: definition.add_optional("name", "string", "name for the representation", default="highres")
else: definition.add_required("name", "string", "name for the representation")

# Name of the model for which to create the representation
model_names = suite.model_names
if len(model_names) == 0: raise RuntimeError("No models found: first run build_model to create a new model")
elif len(model_names) == 1: definition.add_fixed("model_name", "name of the model", model_names[0])
else: definition.add_required("model_name", "string", "name of the model", choices=model_names)

# Whether quality has to be calculated
definition.add_flag("check_dust_grid_quality", "check the quality of the dust grid in various ways", False)

# -----------------------------------------------------------------
