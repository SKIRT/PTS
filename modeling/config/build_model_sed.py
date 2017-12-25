#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.build.suite import ModelSuite
from pts.modeling.core.environment import verify_modeling_cwd
from pts.modeling.config.component import definition

# -----------------------------------------------------------------

modeling_path = verify_modeling_cwd()
suite = ModelSuite.from_modeling_path(modeling_path)

# -----------------------------------------------------------------

definition = definition.copy()

# Add settings
model_names = suite.model_names
if len(model_names) == 0: definition.add_positional_optional("name", "string", "name for the model", default="standard")
else: definition.add_required("name", "string", "name for the model", forbidden=model_names)

# -----------------------------------------------------------------

# Force overwrite
definition.add_flag("overwrite", "overwrite possibly existing model with this name", False)

# -----------------------------------------------------------------
