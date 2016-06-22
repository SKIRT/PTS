#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.photometry Calculate the photometry of the different input images

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.modeling.photometry.photometry import PhotoMeter
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import Configuration

# -----------------------------------------------------------------

# Create the configuration
config = Configuration()

# Read the configuration settings from the provided command-line arguments
config.read()

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = fs.join(fs.cwd(), "log", time.unique_name("log") + ".txt") if config.arguments.report else None

# Determine the log level
level = "DEBUG" if config.arguments.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting photometry ...")

# -----------------------------------------------------------------

# Create a PhotoMeter object
photometer = PhotoMeter(config.get_settings())

# Run the photometer
photometer.run()

# -----------------------------------------------------------------
