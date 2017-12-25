#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.config.build_representation_galaxy import dust_grid_types, default_dust_grid_type
from pts.modeling.build.suite import ModelSuite
from pts.modeling.core.environment import verify_modeling_cwd
from pts.modeling.config.component import definition

# -----------------------------------------------------------------

# Determine the modeling path
modeling_path = verify_modeling_cwd()
suite = ModelSuite.from_modeling_path(modeling_path)

# -----------------------------------------------------------------

default_nrepresentations = 2

# -----------------------------------------------------------------

definition = definition.copy()

# Name of the model for which to create the representation
model_names = suite.model_names
if len(model_names) == 0: raise RuntimeError("No models found: first run build_model to create a new model")
elif len(model_names) == 1: definition.add_fixed("model_name", "name of the model", model_names[0])
else: definition.add_required("model_name", "string", "name of the model", choices=model_names)

# Number of representation to generate
definition.add_optional("nrepresentations", "positive_integer", "number of representations to generate", default_nrepresentations, min_value=2)

# Settings for the dust grid generation
definition.add_section("dg", "settings for the dust grids")
definition.sections["dg"].add_optional("grid_type", "string", "type of dust grid", default_dust_grid_type, choices=dust_grid_types)
definition.sections["dg"].add_optional("scale_range", "real_range", "range of the number of image pixels to take as the minimum scale in the model (can also be a certain fraction of a pixel)", "0.5>10.", convert_default=True)
definition.sections["dg"].add_optional("bintree_level_range", "integer_range", "range of the minimum depth level for binary trees", "6>9", convert_default=True)
definition.sections["dg"].add_optional("octtree_level_range", "integer_range", "range of the minimum depth level for octtrees", "2>3", convert_default=True)
definition.sections["dg"].add_optional("mass_fraction_range", "real_range", "range of the maximum mass fraction in each cell", "0.5e-6>1e-5", convert_default=True)
definition.sections["dg"].add_optional("scale_heights", "real", "number of times to take the dust scale height as the vertical radius of the dust grid", 10.)

# Whether quality has to be calculated
definition.add_flag("check_dust_grid_quality", "check the quality of the dust grid in various ways", True)

# -----------------------------------------------------------------
