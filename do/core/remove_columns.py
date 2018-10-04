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
from pts.core.tools import filesystem as fs
from pts.core.tools import sequences

# -----------------------------------------------------------------

# Create configuration
definition = ConfigurationDefinition(write_config=False)
#definition.add_required("filename", "file_path", "table file")
definition.add_required("columns", "string_list", "names of columns")
definition.add_optional("method", "string", "table reading method", "lines")
config = parse_arguments("remove_columns", definition)

# -----------------------------------------------------------------

# Find table files
filepaths = fs.files_in_cwd(extension="dat", recursive=True)

for filepath in filepaths:

    #print(filepath)
    directory_name = fs.name(fs.directory_of(filepath))
    filename = fs.name(filepath)

    # Get column names
    column_names = fs.get_column_names(filepath)
    if not sequences.contains_any(column_names, config.columns): continue

    print(directory_name, filename, column_names)

    # -----------------------------------------------------------------

    # Load the table
    table = SmartTable.from_file(filepath, method=config.method)

    # -----------------------------------------------------------------

    # Remove the columns
    table.remove_columns(config.columns)

    # -----------------------------------------------------------------

    # Save the table
    table.save()

# -----------------------------------------------------------------
