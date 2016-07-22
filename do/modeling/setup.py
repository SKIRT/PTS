#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.setup Setup a directory for performing the radiative transfer modeling pipeline
#  on a certain galaxy.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import logging
from pts.magic.tools import catalogs
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import ConfigurationDefinition, ConfigurationReader

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

# Add required settings
definition.add_required("name", "string", "the name of the galaxy")

# Get configuration
reader = ConfigurationReader("setup")
config = reader.read(definition)

# -----------------------------------------------------------------

# Determine the log level
level = "DEBUG" if config.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level)
log.start("Starting setup ...")

# -----------------------------------------------------------------

# Inform the user
log.info("Resolving the galaxy name ...")

# Get the NGC name of the galaxy
ngc_name = catalogs.get_ngc_name(config.name)

# Inform the user
log.info("Galaxy NGC ID is '" + ngc_name + "'")

# Determine the path to the new directory
path = fs.join(fs.cwd(), ngc_name)

# Create the directory
fs.create_directory(path)

# Determine the path to the data directory
data_path = fs.join(path, "data")

# Create the data directory
fs.create_directory(data_path)

# -----------------------------------------------------------------
