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

# Determine the modeling path
environment = load_modeling_environment_cwd()
suite = environment.static_model_suite

# -----------------------------------------------------------------

dust_grid_types = ["cartesian", "bintree", "octtree"]
default_dust_grid_type = "bintree"

# -----------------------------------------------------------------

definition = definition.copy()

# Name of the representation
if suite.has_representations: definition.add_required("name", "string", "name for the representation", forbidden=suite.representation_names)
else: definition.add_optional("name", "string", "name for the representation", default="highres")

# Name of the model for which to create the representation
if suite.no_models: raise RuntimeError("No models found: first run build_model to create a new model")
elif suite.has_single_model: definition.add_fixed("model_name", "name of the model", suite.single_model_name)
else: definition.add_required("model_name", "string", "name of the model", choices=suite.model_names)

# Dust grid properties
definition.add_section("dg", "settings for the dust grid")
definition.sections["dg"].add_optional("grid_type", "string", "type of dust grid", default_dust_grid_type, choices=dust_grid_types)
definition.sections["dg"].add_optional("scale", "real", "number of image pixels to take as the minimum scale in the model (can also be a certain fraction of a pixel)", 0.5)
definition.sections["dg"].add_optional("bintree_min_level", "integer", "minimum depth level for binary trees", 9)
definition.sections["dg"].add_optional("octtree_min_level", "integer", "minimum depth level for octrees", 3)
definition.sections["dg"].add_optional("max_mass_fraction", "real", "maximum mass fraction in each cell", 0.5e-6)
definition.sections["dg"].add_optional("scale_heights", "real", "number of times to take the dust scale height as the vertical radius of the dust grid", 10.)

# Whether quality has to be calculated
definition.add_flag("check_dust_grid_quality", "check the quality of the dust grid in various ways", True)

# -----------------------------------------------------------------

# Vertical extent of the total model
definition.add_optional("old_scale_heights", "real", "number of times to take the old stellar scale height as the vertical radius of the model", 15)
definition.add_flag("from_projection", "create the projections for the other orientations from the earth projection (instead of from the deprojection model)", True)
definition.add_optional("radial_factor", "real", "factor with which to multiply the radial extent of the projections", 1.5)

# -----------------------------------------------------------------
