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

# Settings for the wavelength grid generation
definition.add_section("wg", "settings for the wavelength grids")
definition.sections["wg"].add_optional("npoints_range", "integer_range", "the range of the wavelength grid size", "150>500", convert_default=True)
definition.sections["wg"].add_optional("ngrids", "integer", "the number of wavelength grids to generate", 10)
definition.sections["wg"].add_flag("add_emission_lines", "add emission lines to the wavelength grids", False)
definition.sections["wg"].add_optional("min_wavelength", "quantity", "minimum wavelength")
definition.sections["wg"].add_optional("max_wavelength", "quantity", "maximum wavelength")

# Settings for the dust grid generation
definition.add_section("dg", "settings for the dust grids")
definition.sections["dg"].add_optional("grid_type", "string", "the type of dust grid", "bintree", choices=["cartesian", "bintree", "octtree"])
definition.sections["dg"].add_optional("scale_range", "real_range", "the range of the number of image pixels to take as the minimum scale in the model (can also be a certain fraction of a pixel)", "0.5>10.", convert_default=True)
definition.sections["dg"].add_optional("level_range", "integer_range", "the range of the maximum depth level of the tree", "6>9", convert_default=True)
definition.sections["dg"].add_optional("mass_fraction_range", "real_range", "the range of the maximum mass fraction in each cell", "0.5e-6>1e-5", convert_default=True)

# -----------------------------------------------------------------
