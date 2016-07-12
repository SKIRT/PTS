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

# Import the relevant PTS classes and modules
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import Configuration

# -----------------------------------------------------------------

# Create the configuration
config = Configuration("clear")

# Add setting
config.add_required("step", str, "the modeling step for which to clear the output", to_instance=False)

# Read the configuration settings from the provided command-line arguments
config.read()

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = fs.join(fs.cwd(), "log", time.unique_name("log") + ".txt") if config.arguments.report else None

# Determine the log level
level = "DEBUG" if config.arguments.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting clear ...")

# -----------------------------------------------------------------

prep_path = fs.join(config.arguments.path, "prep")
components_path = fs.join(config.arguments.path, "components")
truncated_path = fs.join(config.arguments.path, "truncated")
phot_path = fs.join(config.arguments.path, "phot")
maps_path = fs.join(config.arguments.path, "maps")
fit_path = fs.join(config.arguments.path, "fit")
analysis_path = fs.join(config.arguments.path, "analysis")

if config.arguments.step == "prep": fs.clear_directory(prep_path)
elif config.arguments.step == "components": fs.clear_directory(components_path)
elif config.arguments.step == "truncated": fs.clear_directory(truncated_path)
elif config.arguments.step == "phot": fs.clear_directory(phot_path)
elif config.arguments.step == "maps": fs.clear_directory(maps_path)
elif config.arguments.step == "fit": fs.clear_directory(fit_path)
elif config.arguments.step == "analysis": fs.clear_directory(analysis_path)
else: raise ValueError("Invalid modeling step")

# -----------------------------------------------------------------
