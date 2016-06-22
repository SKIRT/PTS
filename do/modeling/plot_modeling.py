#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.plot_modeling Make plots for a particular radiative transfer modeling step.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import logging, parsing, time
from pts.core.tools import filesystem as fs
from pts.modeling.plotting.plotter import Plotter
from pts.core.basics.configuration import Configuration

# -----------------------------------------------------------------

# Create the configuration
config = Configuration()

# Add required arguments
config.add_required("step", str, "the modeling step for which plots should be made")

# Add optional
config.add_optional("features", "string_list", "the features to be plotted (if not specified this means all features will be plotted")

# Read the configuration settings from the provided command-line arguments
config.read()

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = fs.join(fs.cwd(), "log", time.unique_name("log") + ".txt") if config.arguments.report else None

# Determine the log level
level = "DEBUG" if config.arguments.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting plot_modeling ...")

# -----------------------------------------------------------------

# Create a Plotter instance
plotter = Plotter(config.get_settings())

# Run the plotter
plotter.run()

# -----------------------------------------------------------------
