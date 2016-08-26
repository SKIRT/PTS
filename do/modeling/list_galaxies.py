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
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs
from pts.core.tools import introspection
from pts.dustpedia.core.database import DustPediaDatabase, get_account
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Get configuration
setter = ArgumentConfigurationSetter("list_galaxies")
config = setter.run(definition)

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = fs.join(fs.cwd(), time.unique_name("log") + ".txt") if config.report else None

# Determine the log level
level = "DEBUG" if config.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting list_galaxies ...")

# -----------------------------------------------------------------

# Local table path
local_table_path = fs.join(introspection.pts_dat_dir("modeling"), "s4g", "s4g_p4_table8.dat")

# -----------------------------------------------------------------

def get_galaxy_names_s4g():

    """
    This function ...
    :return:
    """

    names = []

    with open(local_table_path, 'r') as s4g_table:

        for line in s4g_table:

            splitted = line.split()

            if len(splitted) < 2: continue

            # Get the galaxy name and add it to the list
            name = splitted[1]
            names.append(name)

    # Return the list of galaxy names
    return names

# -----------------------------------------------------------------

# Get the account info
username, password = get_account()

# Create the database instance
database = DustPediaDatabase()

# Login with the user and password
#database.login(username, password)

# EARLY TYPE SPIRALS: early-type (Sa–Sab) spiral galaxies

parameters = {"D25": (5., None),
              "Hubble type": ["Sa", "Sab", "Sb"]}

table = database.get_galaxies(parameters)

s4g_names = get_galaxy_names_s4g()

has_s4g_column = []

for i in range(len(table)):

    name = table["Name"][i]

    has_s4g_column.append(name in s4g_names)

table["In S4G"] = has_s4g_column

print(table)

# -----------------------------------------------------------------
