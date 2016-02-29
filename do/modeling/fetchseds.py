#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.fetchseds For a specific galaxy, fetch the SEDs found in various catalogs.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import argparse

# Import the relevant PTS classes and modules
from pts.modeling.photometry.sedfetching import SEDFetcher
from pts.core.tools import logging, time, configuration

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument("name", type=str, help="the name of the galaxy for which to fetch the SEDs")
parser.add_argument("-o", "--output", type=str, help="the name of the output directory")
parser.add_argument("--debug", action="store_true", help="enable debug logging mode")
parser.add_argument("--report", action='store_true', help='write a report file')
parser.add_argument("--config", type=str, help="the name of a configuration file")
parser.add_argument("--settings", type=configuration.from_string, help="settings")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# -- Output --

# If an output directory is given
if arguments.output is not None:

    # Determine the full path to the output directory
    arguments.output_path = os.path.abspath(arguments.output)

    # Create the directory if it does not yet exist
    if not os.path.isdir(arguments.output_path): os.makedirs(arguments.output_path)

# If no output directory is given, place the output in the current working directory
else: arguments.output_path = os.getcwd()

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = os.path.join(arguments.output_path, time.unique_name("sedfetching") + ".txt") if arguments.report else None

# Determine the log level
level = "DEBUG" if arguments.debug else "INFO"

# Initialize the logger
logging.setup_log(level=level, path=logfile_path)
logging.log.info("Starting fetchseds script ...")

# -----------------------------------------------------------------

# Create a SEDFetcher object
fetcher = SEDFetcher.from_arguments(arguments)

# Run the SED fetching
fetcher.run(arguments.name)

# -----------------------------------------------------------------
