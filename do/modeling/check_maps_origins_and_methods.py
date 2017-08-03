#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.check_maps_origins_and_methods Checks ...

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.tools import filesystem as fs
from pts.core.tools.serialization import load_dict
from pts.core.tools import sequences
from pts.core.tools import formatting as fmt

# -----------------------------------------------------------------

# Find methods file
methods_path = fs.join(fs.cwd(), "methods.txt")
if not fs.is_file(methods_path): raise IOError("Methods file not found")
methods = load_dict(methods_path)

# Find origins file
origins_path = fs.join(fs.cwd(), "origins.txt")
if not fs.is_file(origins_path): raise IOError("Origins file not found")
origins = load_dict(origins_path)

# -----------------------------------------------------------------

filenames = fs.files_in_path(fs.cwd(), returns="name", extension="fits")

# -----------------------------------------------------------------

map_names = sequences.union(methods.keys(), origins.keys(), filenames)

# -----------------------------------------------------------------

print("")
print("Presumed number of maps: " + str(len(map_names)))
print("")
print(" - number of map files found: " + str(len(filenames)))
print(" - number of entries in 'origins': " + str(len(origins)))
print(" - number of entries in 'methods': " + str(len(methods)))
print("")

# -----------------------------------------------------------------

# Loop over the maps
for name in map_names:

    # Determine map path
    filepath = fs.join(fs.cwd(), name + ".fits")

    all_found = True

    if not fs.is_file(filepath):
        print(fmt.red + name + ": map not present" + fmt.reset)
        all_found = False

    if not name in origins:
        print(fmt.red + name + ": not present in origins" + fmt.reset)
        all_found = False

    if not name in methods:
        print(fmt.red + name + ": not present in methods" + fmt.reset)
        all_found = False

    # OK
    if all_found: print(fmt.green + name + ": OK" + fmt.reset)

# -----------------------------------------------------------------

print("")

# -----------------------------------------------------------------