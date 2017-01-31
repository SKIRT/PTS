#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.developer.parse Test the parsing of a string into a desired property.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import collections

# Import the relevant PTS classes and modules
from pts.core.tools import parsing
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter

# -----------------------------------------------------------------

definition = ConfigurationDefinition()
definition.add_required("parsing_type", "string", "parsing type")
definition.add_required("string", "string", "string to be parsed")

# -----------------------------------------------------------------

# Get config
setter = ArgumentConfigurationSetter("parse")
config = setter.run(definition)

# -----------------------------------------------------------------

# Get parsing function
parsing_function = getattr(parsing, config.parsing_type)

# Parse
result = parsing_function(config.string)

# Show result
if isinstance(result, collections.Iterable):

    for item in result: print(str(item))

else: print(str(result))

# -----------------------------------------------------------------