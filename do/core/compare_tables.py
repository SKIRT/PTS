#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.compare_tables Compare two tables.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from pts.core.tools import formatting as fmt
from pts.core.tools import tables
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments

# -----------------------------------------------------------------

definition = ConfigurationDefinition(write_config=False)
definition.add_required("path_a", "file_path", "path of the first table")
definition.add_required("path_b", "file_path", "path of the second table")

config = parse_arguments("compare_tables", definition, add_logging=False, add_cwd=False)

# -----------------------------------------------------------------

format_a = "ascii." + fs.get_extension(config.path_a)
format_b = "ascii." + fs.get_extension(config.path_b)

table_a = tables.from_file(config.path_a, format=format_a)
table_b = tables.from_file(config.path_b, format=format_b)

# -----------------------------------------------------------------

print("")

print(fmt.underlined + fmt.green + "Table a" + fmt.reset + " [" + config.path_a + "]:")
print("")
print(" - length: " + str(len(table_a)))
print(" - number of columns: " + str(len(table_a.colnames)))
#print(" - columns: " + stringify.stringify(table_a.colnames)[1])

print("")

print(fmt.underlined + fmt.green + "Table b" + fmt.reset + " [" + config.path_b + "]:")
print("")
print(" - length: " + str(len(table_b)))
print(" - number of columns: " + str(len(table_b.colnames)))
#print(" - columns: " + stringify.stringify(table_b.colnames)[1])
print("")

if table_a.colnames == table_b.colnames: print("Column names are equal")
else:
    print("Column names are not equal")
    differences = list(set(table_a.colnames) - set(table_b.colnames))
    if len(differences) > 0:
        print("Differences:")
        print(differences)
    differences_order = set(table_a.colnames).symmetric_difference(table_b.colnames)
    if len(differences_order) > 0:
        print("Differences in order:")
        print(differences_order)

if len(table_a) == len(table_b): print("Table sizes are equal")
else: print("Table sizes are not equal")

# -----------------------------------------------------------------

for index_a in range(len(table_a)):

    key = table_a[table_a.colnames[0]][index_a]
    index_b = tables.find_index(table_b, key)

    row_a = table_a[index_a]
    row_b = table_b[index_b]

    #print("")
    #print(row_a)
    #print(row_b)
    #print("")

    for name in table_a.colnames:

        try:
            value_a = float(row_a[name])
            value_b = float(row_b[name])
            #print(value_a/value_b)

            if not np.isclose(value_a, value_b): print("DIFFERENCE FOR KEY '" + key + "' and column '" + name + "': " + str(value_a) + " =/= " + str(value_b))

        except: pass

#for index in range(max(len(table_a), len(table_b))):

    #row_a = table_a[index]
    #row_b = table_b[index]

    #if row_a != row_b:
        #print("Row " + str(index) + ":")
        #print(row_a)
        #print(row_b)
        #print("")

    #for name in table_a.colnames:
    #    try:
    #        value_a = float(row_a[name])
    #        value_b = float(row_b[name])
    #        print(value_a/value_b)
    #    except: pass

# -----------------------------------------------------------------
