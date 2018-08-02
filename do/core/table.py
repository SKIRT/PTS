#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.table Show table info.

# -----------------------------------------------------------------

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.tools import filesystem as fs
from pts.core.tools.stringify import tostr
from pts.core.tools import formatting as fmt

# -----------------------------------------------------------------

definition = ConfigurationDefinition()
definition.add_required("filename", "file_path", "filename or filepath of the data")
definition.add_positional_optional("nlines", "positive_integer", "number of first and last lines to show", 5)
config = parse_arguments("table", definition)

# -----------------------------------------------------------------

print("")

# -----------------------------------------------------------------

print(fmt.bold + "HEADER:" + fmt.reset_bold)
print("")

header = fs.get_header_lines(config.filename)
nheaderlines = len(header)
for line in header: print(line)

# -----------------------------------------------------------------

print("")

# -----------------------------------------------------------------

colnames = fs.get_column_names(config.filename)
print(fmt.bold + "COLUMN NAMES: " + fmt.reset_bold + tostr(colnames))

# -----------------------------------------------------------------

print("")

# -----------------------------------------------------------------

nrows = fs.get_nlines(config.filename) - nheaderlines
print(fmt.bold + "NUMBER OF ROWS: " + fmt.reset_bold + str(nrows))

# -----------------------------------------------------------------

print("")

# -----------------------------------------------------------------

print(fmt.bold + "FIRST " + str(config.nlines) + " LINES:" + fmt.reset_bold)
print("")

first = fs.get_first_lines(config.filename, config.nlines + nheaderlines)
for index, line in enumerate(first):
    if index < nheaderlines: continue
    print(line)

# -----------------------------------------------------------------

print("")

# -----------------------------------------------------------------

print(fmt.bold + "LAST " + str(config.nlines) + " LINES:" + fmt.reset_bold)
print("")

last = fs.get_last_lines(config.filename, config.nlines)
for line in last: print(line)

# -----------------------------------------------------------------

print("")

# -----------------------------------------------------------------
