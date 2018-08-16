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
from pts.core.basics.log import log
from pts.core.tools import parsing
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments

# -----------------------------------------------------------------

definition = ConfigurationDefinition()
definition.add_required("parsing_type", "string", "parsing type")
definition.add_required("string", "string", "string to be parsed")

# -----------------------------------------------------------------

# Get config
config = parse_arguments("parse", definition)

# -----------------------------------------------------------------

# Get parsing function
parsing_function = getattr(parsing, config.parsing_type)

# Parse
try:

    result = parsing_function(config.string)

    # Show result
    if isinstance(result, collections.Iterable):
        for item in result: print(str(item))
    else: print(str(result))

except ValueError as e:
    log.error("The string could not be parsed into this property")
    log.error(e.message)

# -----------------------------------------------------------------
