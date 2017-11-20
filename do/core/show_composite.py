#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.show_composite Show a property composite

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.basics.composite import load_composite
from pts.core.basics.table import SmartTable

# -----------------------------------------------------------------

definition = ConfigurationDefinition(write_config=False)
definition.add_required("filename", "file_path", "path of the property composite file")
definition.add_flag("latex", "print as latex")

config = parse_arguments("show_composite", definition, add_logging=False, add_cwd=False)

# -----------------------------------------------------------------

# Load
composite = load_composite(config.filename)

# -----------------------------------------------------------------

if config.latex:

    table = SmartTable.from_composite(composite)
    table.print_latex()

else: print(composite)

# -----------------------------------------------------------------
