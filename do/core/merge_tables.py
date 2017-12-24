#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.merge_tables Merge tables into one.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.basics.table import SmartTable, merge_tables
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()
definition.add_required("files", "filepath_list", "list of the table file paths")
definition.add_positional_optional("column", "string", "name of column to use as reference for the merge")
definition.add_flag("show", "show the merged table", False)
definition.add_flag("write", "write the merged table", True)
definition.add_flag("differences", "add absolute differences", False)
definition.add_flag("rel_differences", "add relative differences", False)
definition.add_flag("percentual", "put relative differences in percentages", False)
definition.add_optional("columns", "string_list", "only use these columns (other than the reference column)")
definition.add_optional("not_columns", "string_list", "don't use these columns")

# Parse
config = parse_arguments("merge_tables", definition)

# -----------------------------------------------------------------

# Load the tables
tables = []
for filepath in config.files:
    table = SmartTable.from_file(filepath)
    tables.append(table)

# -----------------------------------------------------------------

filenames = [fs.strip_extension(fs.name(filepath)) for filepath in config.files]

# -----------------------------------------------------------------

table = merge_tables(*tables, column_name=config.column, labels=filenames, differences=config.differences,
                     rel_differences=config.rel_differences, percentual=config.percentual, columns=config.columns, not_columns=config.not_columns)

# -----------------------------------------------------------------

# Show the table
if config.show: print(table)

# Write the table
if config.write:

    # Determine path
    path = "merged.dat"

    # Save
    table.saveto(path)

# -----------------------------------------------------------------
