#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.dustpedia.core.database import DustPediaDatabase, get_account
from pts.dustpedia.core.sample import DustPediaSample
from pts.core.tools import filesystem as fs
from pts.core.tools.stringify import tostr
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments

# -----------------------------------------------------------------

definition = ConfigurationDefinition()
definition.add_required("galaxy_name", "string", "galaxy name")

# Set
config = parse_arguments("get_mbb_parameters", definition)

# -----------------------------------------------------------------

path = fs.cwd()

# -----------------------------------------------------------------

# Create the DustPedia sample
sample = DustPediaSample()
galaxy_name = sample.get_name(config.galaxy_name)

# Create the database
database = DustPediaDatabase()

# Login
username, password = get_account()
database.login(username, password)

# Get the parameters
parameters = database.get_dust_black_body_parameters(galaxy_name)

dust_mass = parameters[0]
dust_mass_error = parameters[1]

temperature = parameters[2]
temperature_error = parameters[3]

luminosity = parameters[4]

print("")
print(galaxy_name)
print(" - Mass: " + tostr(dust_mass, scientific=True, ndigits=3))
print(" - Temperature: " + tostr(temperature, scientific=True, ndigits=3))
print(" - Luminosity: " + tostr(luminosity, scientific=True, ndigits=3))
print("")

# -----------------------------------------------------------------
