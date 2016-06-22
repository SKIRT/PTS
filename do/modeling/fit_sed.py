#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.fit_sed Determine the best fit with the observed SED

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.modeling.fitting.sedfitting import SEDFitter
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import Configuration

# -----------------------------------------------------------------

# Create the configuration
config = Configuration()

# Add optional arguments
config.add_flag("visualise", "make visualisations")

# Read the configuration settings from the provided command-line arguments
config.read()

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = fs.join(fs.cwd(), "log", time.unique_name("log") + ".txt") if config.arguments.report else None

# Determine the log level
level = "DEBUG" if config.arguments.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting fit_sed ...")

# -----------------------------------------------------------------

# Create a SEDFitter object
fitter = SEDFitter(config.get_settings())

# Run the fitter
fitter.run()

# -----------------------------------------------------------------
