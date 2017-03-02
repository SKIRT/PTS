#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.developer.configure_skirt_item Test the configuration of a SKIRT simulation item.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.prep.smile import SKIRTSmileSchema
from pts.core.basics.configuration import ConfigurationDefinition, print_mapping
from pts.core.basics.configuration import ArgumentConfigurationSetter
from pts.core.tools import stringify
from pts.core.tools import formatting as fmt

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()
definition.add_positional_optional("name", "string", "name of the configurable item")

# Parse
setter = ArgumentConfigurationSetter("configure_skirt_item")
config = setter.run(definition)

# -----------------------------------------------------------------

# Create the SKIRT smile schema
smile = SKIRTSmileSchema()

# Get the configuration parameters interactively
parameters, children = smile.prompt_parameters_for_type(config.name)

print("")
print(fmt.red + fmt.underlined + "Parameters" + fmt.reset)

# Print mapping
print_mapping(parameters)

# Show simulation items
print(fmt.blue + "Simulation items: " + fmt.reset + stringify.stringify(children.keys())[1].replace(",", ", "))
print("")

# Print children
for name in children:

    # Get parameters of child (and children)
    parameters, child_children = children[name]

    if len(parameters) == 0:

        print("    " + fmt.red + fmt.underlined + name + ": no parameters" + fmt.reset)
        print("")

    else:

        print("    " + fmt.red + fmt.underlined + name + " parameters" + fmt.reset)

        # Print mapping
        print_mapping(parameters, indent="    ")

        # Show simulation items

#print("")

# -----------------------------------------------------------------
