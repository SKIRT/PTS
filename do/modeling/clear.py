#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.clear Clear the output of one of the radiative transfer modeling steps.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import argparse

# Import the relevant PTS classes and modules
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()

# Basic options
parser.add_argument("step", type=str, help="the preparation step for which to clear the output")

# Logging
parser.add_argument("--debug", action="store_true", help="enable debug logging mode")
parser.add_argument("--report", action='store_true', help='write a report file')

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Set the modeling path and the log path
arguments.path = fs.cwd()
log_path = fs.join(arguments.path, "log")

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = fs.join(log_path, time.unique_name("log") + ".txt") if arguments.report else None

# Determine the log level
level = "DEBUG" if arguments.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting clear ...")

# -----------------------------------------------------------------

prep_path = fs.join(arguments.path, "prep")
components_path = fs.join(arguments.path, "components")
truncated_path = fs.join(arguments.path, "truncated")
phot_path = fs.join(arguments.path, "phot")
maps_path = fs.join(arguments.path, "maps")
fit_path = fs.join(arguments.path, "fit")
analysis_path = fs.join(arguments.path, "analysis")

if arguments.step == "prep": fs.clear_directory(prep_path)
elif arguments.step == "components": fs.clear_directory(components_path)
elif arguments.step == "truncated": fs.clear_directory(truncated_path)
elif arguments.step == "phot": fs.clear_directory(phot_path)
elif arguments.step == "maps": fs.clear_directory(maps_path)
elif arguments.step == "fit": fs.clear_directory(fit_path)
elif arguments.step == "analysis": fs.clear_directory(analysis_path)
else: raise ValueError("Invalid modeling step")

# -----------------------------------------------------------------
