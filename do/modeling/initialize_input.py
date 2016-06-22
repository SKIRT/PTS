#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.initialize_input Initialize the input directory for the fitting procedure.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.modeling.fitting.initialization import InputInitializer
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import Configuration

# -----------------------------------------------------------------

# Create the configuration
config = Configuration()

# Add optional arguments
config.add_optional("nlambda", "make visualisations")
config.add_optional("lambda_minmax", "float_tuple", "the minimum and maximum wavelength of the total grid")
config.add_optional("lambda_minmax_zoom", "float_tuple", "the minimum and maximum wavelength of the zoomed-in grid")
config.add_optional("packages", float, "the number of photon packages per wavelength")
config.add_flag("selfabsorption", "enable dust self-absorption")

# Read the configuration settings from the provided command-line arguments
config.read()

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = fs.join(fs.cwd(), "log", time.unique_name("initialization") + ".txt") if config.arguments.report else None

# Determine the log level
level = "DEBUG" if config.arguments.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting initialize_input ...")

# -----------------------------------------------------------------

# Create a InputInitializer object
initializer = InputInitializer(config.get_settings())

# Run the input initialization
initializer.run()

# -----------------------------------------------------------------
