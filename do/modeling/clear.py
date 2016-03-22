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
from pts.core.tools import logging, time, filesystem

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

modeling_path = filesystem.cwd()

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = filesystem.join(modeling_path, time.unique_name("log") + ".txt") if arguments.report else None

# Determine the log level
level = "DEBUG" if arguments.debug else "INFO"

# Initialize the logger
logging.setup_log(level=level, path=logfile_path)
logging.log.info("Starting clear ...")

# -----------------------------------------------------------------

prep_path = filesystem.join(modeling_path, "prep")
components_path = filesystem.join(modeling_path, "components")
truncated_path = filesystem.join(modeling_path, "truncated")
phot_path = filesystem.join(modeling_path, "phot")
maps_path = filesystem.join(modeling_path, "maps")
fit_path = filesystem.join(modeling_path, "fit")
analysis_path = filesystem.join(modeling_path, "analysis")

if arguments.step == "prep": filesystem.clear_directory(prep_path)
elif arguments.step == "components": filesystem.clear_directory(components_path)
elif arguments.steps == "truncated": filesystem.clear_directory(truncated_path)
elif arguments.steps == "phot": filesystem.clear_directory(phot_path)
elif arguments.steps == "maps": filesystem.clear_directory(maps_path)
elif arguments.steps == "fit": filesystem.clear_directory(fit_path)
elif arguments.steps == "analysis": filesystem.clear_directory(analysis_path)
else: raise ValueError("Invalid modeling step")

# -----------------------------------------------------------------
