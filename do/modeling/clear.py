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
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add setting
definition.add_required("step", "string", "the modeling step for which to clear the output")

# Get the configuration
setter = ArgumentConfigurationSetter("clear")
config = setter.run(definition)

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = fs.join(fs.cwd(), "log", time.unique_name("log") + ".txt") if config.arguments.report else None

# Determine the log level
level = "DEBUG" if config.arguments.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting clear ...")

# -----------------------------------------------------------------

prep_path = fs.join(config.path, "prep")
components_path = fs.join(config.path, "components")
truncated_path = fs.join(config.path, "truncated")
phot_path = fs.join(config.path, "phot")
maps_path = fs.join(config.path, "maps")
fit_path = fs.join(config.path, "fit")
analysis_path = fs.join(config.path, "analysis")

if config.step == "prep": fs.clear_directory(prep_path)
elif config.step == "components": fs.clear_directory(components_path)
elif config.step == "truncated": fs.clear_directory(truncated_path)
elif config.step == "phot": fs.clear_directory(phot_path)
elif config.step == "maps": fs.clear_directory(maps_path)
elif config.step == "fit": fs.clear_directory(fit_path)
elif config.step == "analysis": fs.clear_directory(analysis_path)
else: raise ValueError("Invalid modeling step")

# -----------------------------------------------------------------
