#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.developer.skirt_units List the SKIRT unit systems.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.prep.smile import SKIRTSmileSchema
from pts.core.tools import formatting as fmt
from pts.core.units.stringify import stringify_unit
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.basics.configuration import ArgumentConfigurationSetter

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()
definition.add_positional_optional("match", "string", "only show unit systems with names that contain this string")

# Parse
setter = ArgumentConfigurationSetter("skirt_units")
config = setter.run(definition)

# -----------------------------------------------------------------

# Create the SKIRT smile schema
smile = SKIRTSmileSchema()

print("")

# Loop over
for unit_system in smile.unit_systems:

    # Get default units
    units = smile.units_for_unit_system(unit_system)

    print(fmt.green + fmt.underlined + unit_system + fmt.reset)
    print("")

    for quantity_name in units:

        print(" - " + fmt.blue + quantity_name + fmt.reset + ": " + stringify_unit(units[quantity_name])[1])

    print("")

# -----------------------------------------------------------------
