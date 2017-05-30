#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.modeling.build.component import get_model_names, get_representation_names
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

# Determine the modeling path
modeling_path = fs.cwd()

# -----------------------------------------------------------------

dust_grid_types = ["cartesian", "bintree", "octtree"]
default_dust_grid_type = "bintree"

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# Name of the representation
representation_names = get_representation_names(modeling_path)
if len(representation_names) == 0: definition.add_optional("name", "string", "name for the representation", default="highres")
else: definition.add_required("name", "string", "name for the representation")

# Name of the model for which to create the representation
model_names = get_model_names(modeling_path)
if len(model_names) == 0: raise RuntimeError("No models found: first run build_model to create a new model")
elif len(model_names) == 1: definition.add_fixed("model_name", "name of the model", model_names[0])
else: definition.add_required("model_name", "string", "name of the model", choices=model_names)

# Dust grid properties
definition.add_section("dg", "settings for the dust grid")
definition.sections["dg"].add_optional("grid_type", "string", "type of dust grid", default_dust_grid_type, choices=dust_grid_types)
definition.sections["dg"].add_optional("scale", "real", "number of image pixels to take as the minimum scale in the model (can also be a certain fraction of a pixel)", 0.5)
definition.sections["dg"].add_optional("max_level", "integer", "maximum depth level of the tree", 9)
definition.sections["dg"].add_optional("mass_fraction", "real", "maximum mass fraction in each cell", 0.5e-6)

# -----------------------------------------------------------------
