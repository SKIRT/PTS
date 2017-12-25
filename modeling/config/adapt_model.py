#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.modeling.decomposition.decomposition import scalelength_scaleheight_ratios, degeyter_ratio
from pts.modeling.config.component import definition

# -----------------------------------------------------------------

# Load environment and model suite
environment = load_modeling_environment_cwd()
suite = environment.static_model_suite

# -----------------------------------------------------------------

dust_and_stellar = ["dust", "stellar"]

# -----------------------------------------------------------------

definition = definition.copy()

# No models?
if not suite.has_models: raise ValueError("There are currently no models")

# The model name
definition.add_required("name", "string", "name of the model to be adapted", choices=suite.model_names)

# No components (only representations)
definition.add_flag("no_components", "don't adapt model components, only model representation(s)", False)

# Dust or stellar
definition.add_positional_optional("dust_or_stellar", "string_list", "adapt dust or stellar component(s)", default=dust_and_stellar, choices=dust_and_stellar)

# Name of dust or stellar component to adapt
definition.add_positional_optional("component_name", "string", "name of the dust/stellar component to adapt (only when either 'stellar' or 'dust' is selected)")

# Property name
definition.add_positional_optional("matching", "string", "only adapt properties with a name matching this string")

# -----------------------------------------------------------------

# Also adapt representations
if suite.has_representations: definition.add_optional("representations", "string_list", "adapt properties of these representations", choices=suite.representation_names)
else: definition.add_fixed("representations", "adapt properties of these representations", None)

# -----------------------------------------------------------------

# Save
definition.add_flag("save", "save adapted properties", True)

# -----------------------------------------------------------------

# Select certain properties
definition.add_optional("contains", "string", "only adapt properties containing this string in their name")
definition.add_optional("not_contains", "string", "don't adapt properties containing this string in their name")
definition.add_optional("exact_name", "string", "only adapt properties with this exact string as their name")
definition.add_optional("exact_not_name", "string", "don't adapt properties with this exact string as their name")
definition.add_optional("startswith", "string", "only adapt properties whose name starts with this string")
definition.add_optional("endswith", "string", "only adapt properties whose name starts with this string")

# -----------------------------------------------------------------

# Show after adapting
definition.add_flag("show", "show the components after the model is built", True)

# -----------------------------------------------------------------

definition.add_flag("use_defaults", "use default parameter values", False)

# -----------------------------------------------------------------

definition.add_optional("metallicity", "positive_real", "metallicity for templates") # no default or it is used!

# -----------------------------------------------------------------

# Star formation (Mappings)
definition.add_optional("default_ionizing_compactness", "real", "compactness", 6.)
definition.add_optional("default_ionizing_pressure", "quantity", "pressure", "1e12 K/m3", convert_default=True)
definition.add_optional("default_covering_factor", "real", "covering factor", 0.2)

# -----------------------------------------------------------------

# Scale heights
definition.add_optional("scalelength_to_scaleheight", "real", "ratio of scalelength to scaleheight", default=degeyter_ratio, suggestions=scalelength_scaleheight_ratios)
definition.add_optional("young_scaleheight_ratio", "real", "ratio of the young stellar scaleheight to the old stellar scaleheight", 0.5)
definition.add_optional("ionizing_scaleheight_ratio", "real", "ratio of the ionizing scaleheight to the old stellar scaleheight", 0.25)
definition.add_optional("dust_scaleheight_ratio", "real", "ratio of the dust scaleheight to the old stellar scaleheight", 0.5)

# -----------------------------------------------------------------
