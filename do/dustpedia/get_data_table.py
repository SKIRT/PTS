#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.dustpedia.core.database import get_account
from pts.dustpedia.core.properties import DustPediaProperties

# -----------------------------------------------------------------

definition = ConfigurationDefinition()
definition.add_required("galaxy_name", "string", "galaxy name")
config = parse_arguments("get_data_table", definition)

# -----------------------------------------------------------------

# Get DustPedia username and password
username, password = get_account()

# -----------------------------------------------------------------

properties = DustPediaProperties()

# Get the table
table = properties.create_table(config.galaxy_name, username, password)

table.print_latex()

# -----------------------------------------------------------------
