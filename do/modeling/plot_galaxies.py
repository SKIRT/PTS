#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.plot_galaxies Plot the positions of the galaxies in the DustPedia database.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.core.tools import introspection
from pts.dustpedia.core.database import DustPediaDatabase, get_account
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

# Get configuration
setter = parse_arguments("plot_galaxies", definition)

# -----------------------------------------------------------------

# Local table path
local_table_path = fs.join(introspection.pts_dat_dir("modeling"), "s4g", "s4g_p4_table8.dat")

# -----------------------------------------------------------------

# Get the account info
username, password = get_account()

# Create the database instance
database = DustPediaDatabase()

# Login with the user and password
#database.login(username, password)

table = database.get_galaxies(parameters)

s4g_names = get_galaxy_names_s4g()

has_s4g_column = []

for i in range(len(table)):

    name = table["Name"][i]

    has_s4g_column.append(name in s4g_names)

table["In S4G"] = has_s4g_column

print(table)

# -----------------------------------------------------------------
