#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.remove_columns Remove column(s) from a table.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.table import SmartTable
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments

# -----------------------------------------------------------------

# Create configuration
definition = ConfigurationDefinition(write_config=False)
definition.add_required("filename", "file_path", "table file")
definition.add_required("columns", "string_list", "names of columns")
definition.add_optional("method", "string", "table reading method", "lines")
config = parse_arguments("remove_columns", definition)

# -----------------------------------------------------------------

# Load the table
table = SmartTable.from_file(config.filename, method=config.method)

# -----------------------------------------------------------------

# Remove the columns
table.remove_columns(config.columns)

# -----------------------------------------------------------------

# Save the table
table.save()

# -----------------------------------------------------------------
