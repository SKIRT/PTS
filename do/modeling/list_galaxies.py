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

# Import standard modules
import argparse

# Import the relevant PTS classes and modules
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs
from pts.magic.misc.dustpedia import DustPediaDatabase, get_account

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()

# Logging
parser.add_argument("--debug", action="store_true", help="enable debug logging mode")
parser.add_argument("--report", action='store_true', help='write a report file')

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = fs.join(arguments.path, time.unique_name("log") + ".txt") if arguments.report else None

# Determine the log level
level = "DEBUG" if arguments.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting list_galaxies ...")

# -----------------------------------------------------------------

# Get the account info
username, password = get_account()

# Create the database instance
database = DustPediaDatabase()

# Login with the user and password
database.login(username, password)

#

# -----------------------------------------------------------------
