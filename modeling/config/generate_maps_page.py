#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.config.generate_page import definition
from pts.modeling.build.suite import ModelSuite
from pts.modeling.core.environment import verify_modeling_cwd

# -----------------------------------------------------------------

# Set the modeling path
modeling_path = verify_modeling_cwd()
suite = ModelSuite.from_modeling_path(modeling_path)

# -----------------------------------------------------------------

# Fitting run setting
if suite.no_models: raise RuntimeError("No models are present (yet)")
elif suite.has_single_model: definition.add_fixed("model_name", "name of the model", suite.single_model_name)
else: definition.add_required("model_name", "string", "name of the model", choices=suite.model_names)

# -----------------------------------------------------------------
