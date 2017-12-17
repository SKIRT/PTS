#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.show_dictionary Show a saved dictionary.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.tools.serialization import load_dict
from pts.core.basics.table import SmartTable
from pts.core.tools import formatting as fmt
from pts.core.tools.stringify import tostr

# -----------------------------------------------------------------

definition = ConfigurationDefinition(write_config=False)
definition.add_required("filename", "file_path", "path of the dictionary file")
definition.add_flag("latex", "print as latex")

config = parse_arguments("show_dictionary", definition, add_logging=False, add_cwd=False)

# -----------------------------------------------------------------

# Load
dictionary = load_dict(config.filename)

# -----------------------------------------------------------------

if config.latex:

    table = SmartTable.from_dictionary(dictionary)
    table.print_latex()

else:
    bullet = "-"
    for label in dictionary:
        line = " " + bullet + " " + fmt.bold + label + fmt.reset + ": " + tostr(dict[label])

# -----------------------------------------------------------------
