#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.developer.skirt_items List each configurable SKIRT item, based on the SKIRT smile schema.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import StringIO

# Import the relevant PTS classes and modules
from pts.core.prep.smile import SKIRTSmileSchema
from pts.core.tools import formatting as fmt
from pts.core.tools import stringify
from pts.core.basics.configuration import ConfigurationDefinition, write_definition
from pts.core.basics.configuration import ArgumentConfigurationSetter

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()
definition.add_positional_optional("match", "string", "only show types with names that contain this string")
definition.add_flag("definitions", "format the properties for each item as a configuration definition")

# Parse
setter = ArgumentConfigurationSetter("skirt_items")
config = setter.run(definition)

# -----------------------------------------------------------------

# Create the SKIRT smile schema
smile = SKIRTSmileSchema()

print("")

# Loop over the concrete types
for name in smile.concrete_types:

    if config.match is not None and config.match not in name: continue

    print(fmt.green + fmt.underlined + name + fmt.reset)
    print("")

    description = smile.concrete_types[name]
    print(fmt.darkgray + stringify.stringify_string_fancy(description, lines_prefix=" ")[1] + fmt.reset)
    print("")

    if config.definitions:

        # Create definition
        definition, simulation_items = smile.definition_for_type(name)

        # Create output string
        output = StringIO.StringIO()

        # Write definition to string buffer
        write_definition(definition, output, indent=" ")

        # Print the definition
        contents = output.getvalue()
        line_description = None
        for line in contents.split("\n"):
            if line.strip().startswith("#"):
                line_description = line.split("# ")[1]
                continue
            else:
                if line_description is not None: print(line + "   " + fmt.darkgray + line_description + fmt.reset)
                else: print(line)
                line_description = None

        # Close string buffer
        output.close()

    else:

        # Get properties
        properties = smile.properties_for_type(name)

        # Show the properties
        for prop_name in properties:

            description = properties[prop_name].description
            ptype = properties[prop_name].ptype
            min_value = properties[prop_name].min
            max_value = properties[prop_name].max
            default = properties[prop_name].default
            choices = properties[prop_name].choices
            item = properties[prop_name].item

            string = fmt.blue + prop_name + fmt.reset

            string += " [type: " + ptype + "]"
            if min_value is not None: string += " [min: " + stringify.stringify_not_list(min_value)[1] + "]"
            if max_value is not None: string += " [max: " + stringify.stringify_not_list(max_value)[1] + "]"
            if default is not None: string += " [default: " + stringify.stringify_not_list(default)[1] + "]"
            if choices is not None: string += " [choices: " + stringify.stringify(choices.keys())[1] + "]"

            string += " " + fmt.darkgray + stringify.stringify_string_fancy(description, lines_prefix="    ")[1] + fmt.reset

            if item: string += fmt.red + " [SIMULATION ITEM] " + fmt.reset

            print(" -  " + string)
            print("")

    #print("")

# -----------------------------------------------------------------
