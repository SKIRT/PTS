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
from pts.modeling.build.models.galaxy import show_component
from pts.core.tools import formatting as fmt
from pts.modeling.build.models.stars import bulge_component_name, old_component_name, young_component_name, ionizing_component_name
from pts.modeling.build.models.dust import disk_component_name

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

# Load the model definition
model = suite.get_model_definition(config.model_name)

# -----------------------------------------------------------------

print("")
print(fmt.green + fmt.underlined + "STELLAR COMPONENTS" + fmt.reset)
print("")

# -----------------------------------------------------------------

print("  " + fmt.magenta + "OLD STELLAR BULGE:" + fmt.reset)
print("")

# Show
bulge_path = model.get_stellar_component_path(bulge_component_name)
show_component(bulge_path, line_prefix="    ")
print("")

print("  " + fmt.magenta + "OLD STELLAR DISK: " + fmt.reset)
print("")

# Show
old_path = model.get_stellar_component_path(old_component_name)
show_component(old_path, line_prefix="    ")
print("")

print("  " + fmt.magenta + "YOUNG STELLAR DISK:" + fmt.reset)
print("")

# Show
young_path = model.get_stellar_component_path(young_component_name)
show_component(young_path, line_prefix="    ")
print("")

print("  " + fmt.magenta + "IONIZING STELLAR DISK: " + fmt.reset)
print("")

# Show
ionizing_path = model.get_stellar_component_path(ionizing_component_name)
show_component(ionizing_path, line_prefix="    ")
print("")

# Loop over the extra stellar components
for name in model.additional_stellar_names:

    # Show the name
    print("  " + fmt.magenta + name.upper() + ":" + fmt.reset)
    print("")

    path = model.get_stellar_component_path(name)

    # Show the component
    show_component(path, line_prefix="    ")
    print("")

# -----------------------------------------------------------------

print("")
print(fmt.green + fmt.underlined + "DUST COMPONENTS" + fmt.reset)
print("")

# -----------------------------------------------------------------

print("  " + fmt.magenta + "DUST DISK:" + fmt.reset)
print("")

# Show component
dust_path = model.get_dust_component_path(disk_component_name)
show_component(dust_path, line_prefix="    ")
print("")

# Loop over the extra dust components
for name in model.additional_dust_names:

    # Show the name
    print("  " + fmt.magenta + name.upper() + ": " + fmt.reset)
    print("")

    path = model.get_dust_component_path(name)

    # Show
    show_component(path, line_prefix="    ")
    print("")

# -----------------------------------------------------------------
