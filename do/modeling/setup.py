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

# Import standard modules
import argparse

# Import the relevant PTS classes and modules
from pts.core.tools import logging, time, parsing
from pts.magic.tools import catalogs
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()

# The name of the galaxy
parser.add_argument("name", type=str, help="the name of the galaxy")

# Logging
parser.add_argument("--debug", action="store_true", help="enable debug logging mode")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Determine the log level
level = "DEBUG" if arguments.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level)
log.start("Starting setup ...")

# -----------------------------------------------------------------

# Inform the user
log.info("Resolving the galaxy name ...")

# Get the NGC name of the galaxy
ngc_name = catalogs.get_ngc_name(arguments.name)

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
