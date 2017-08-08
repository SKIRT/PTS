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
from pts.core.basics.log import log
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments

# -----------------------------------------------------------------

# Create configuration
definition = ConfigurationDefinition()
definition.add_required("galaxy_name", "string", "name of the galaxy (will be resolved)")
config = parse_arguments("get_headers", definition)

# -----------------------------------------------------------------

# Create the DustPedia sample
sample = DustPediaSample()
galaxy_name = sample.get_name(config.galaxy_name)

# Create the database
database = DustPediaDatabase()

# Login
username, password = get_account()
database.login(username, password)

# Loop over the images for M81
filters = database.get_image_names_and_filters(galaxy_name)
for name in filters:

    # Get the filter
    fltr = filters[name]

    # Inform the user
    log.info("Getting the header for the '" + str(fltr) + "' filter ...")

    # Get the header
    header = database.get_header(galaxy_name, name)

    # Determine path
    header_path = fs.join(path, str(fltr) + ".txt")

    # Save the header
    header.totextfile(header_path)

# -----------------------------------------------------------------
