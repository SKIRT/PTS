#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.show_model Show the composition of a certain model.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.modeling.misc.examination import ModelExamination

# -----------------------------------------------------------------

# Determine the modeling path
environment = load_modeling_environment_cwd()
suite = environment.static_model_suite

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Name of the model for which to create the representation
if suite.no_models: raise RuntimeError("No models found: first run build_model to create a new model")
elif suite.has_single_model: definition.add_fixed("model_name", "name of the model", suite.single_model_name)
else: definition.add_required("model_name", "string", "name of the model", choices=suite.model_names)

# Get configuration
config = parse_arguments("show_model", definition)

# -----------------------------------------------------------------

# Load the model
model = suite.get_model(config.model_name)

# -----------------------------------------------------------------

# Examination
examination = ModelExamination()

# Set flags
examination.config.show_components = True

# Run
examination.run(model=model)

# -----------------------------------------------------------------
