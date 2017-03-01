#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.developer.skirt_quantities List the SKIRT quantities.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.prep.smile import SKIRTSmileSchema, skirt_quantities_to_pts_quantities
from pts.core.tools import formatting as fmt
from pts.core.tools import stringify
from pts.core.basics.configuration import ConfigurationDefinition, write_definition
from pts.core.basics.configuration import ArgumentConfigurationSetter

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()
definition.add_positional_optional("match", "string", "only show quantities with names that contain this string")

# Parse
setter = ArgumentConfigurationSetter("skirt_quantities")
config = setter.run(definition)

# -----------------------------------------------------------------

# Create the SKIRT smile schema
smile = SKIRTSmileSchema()

print("")
for quantity in smile.quantities:

    units = smile.units_for_quantity(quantity)
    print(fmt.green + fmt.underlined + quantity + fmt.reset + ": " + stringify.stringify(units)[1] + " [" + skirt_quantities_to_pts_quantities[quantity] + "]")
    print("")

#print("")

# -----------------------------------------------------------------
