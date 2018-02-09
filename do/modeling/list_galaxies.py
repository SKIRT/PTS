#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.list_galaxies List the galaxies in the DustPedia database that are eligible for the
#  radiative transfer modeling.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.basics.table import SmartTable
from pts.core.tools import tables

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Required
definition.add_required("filename", "file_path", "data table filename")

# Get configuration
config = parse_arguments("list_galaxies", definition)

# -----------------------------------------------------------------

# Determine full path
filepath = fs.absolute_or_in_cwd(config.filename)

# Load the table
table = SmartTable.from_file(filepath)
galaxy_names = list(table["name"])

# -----------------------------------------------------------------

#print(table.colnames)

# -----------------------------------------------------------------

table.sort("d25")

#index = tables.find_index(table, "NGC3031", "name")
#d25_m81 = table["d25"][index]
#print(d25_m81)

# -----------------------------------------------------------------

#print(table["d25"])

#print(table)

# -----------------------------------------------------------------

names = []

index = len(table)
while len(names) < 10:

    index = index - 1
    #print(index)
    name = table["name"][index]
    has_s4g = table["has_s4g"][index]
    has_galex = table["has_galex"][index]
    has_pacs = table["has_pacs"][index]
    has_spire = table["has_spire"][index]

    if not has_s4g: continue
    if not has_galex: continue
    #if not has_pacs: continue
    #if not has_spire: continue

    inclination = table["inclination"][index]
    if inclination > 65: continue

    #print(name, has_s4g, has_galex, has_pacs, has_spire, inclination)

    common_name = table["common_name"][index]

    print(name, common_name)

    names.append(name)

# -----------------------------------------------------------------
