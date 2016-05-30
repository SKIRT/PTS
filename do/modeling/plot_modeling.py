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

# Import standard modules
import argparse

# Import the relevant PTS classes and modules
from pts.core.tools import logging, parsing, time
from pts.core.tools import filesystem as fs
from pts.modeling.plotting.plotter import Plotter

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()

# The modeling step
parser.add_argument("step", type=str, help="the modeling step for which plots should be made")

# Logging options
parser.add_argument("--debug", action="store_true", help="enable debug logging mode")
parser.add_argument("--report", action='store_true', help='write a report file')
parser.add_argument("--config", type=str, help="the name of a configuration file")

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
log.start("Starting plot_modeling ...")

# -----------------------------------------------------------------

# Create a Plotter instance
plotter = Plotter.from_arguments(arguments)

# Run the plotter
plotter.run()

# -----------------------------------------------------------------
