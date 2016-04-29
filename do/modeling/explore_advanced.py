#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.explore_advanced Explore the parameter space for the radiative transfer modeling.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import argparse

# Import the relevant PTS classes and modules
from pts.modeling.fitting.advancedparameterexplorer import AdvancedParameterExplorer
from pts.core.tools import logging, time, parsing, filesystem

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()

# Basic
parser.add_argument("path", type=str, nargs='?', help="the modeling path")

# The number of simulations to launch in the batch
parser.add_argument("--simulations", type=int, help="the number of simulations to launch in the batch")

# The remote host(s)
parser.add_argument("--remotes", type=parsing.string_list, help="the remote hosts on which to launch the simulations")

# The range for the different fit parameters
parser.add_argument("--young", type=parsing.int_tuple, help="the range of the FUV luminosity of the young stellar population")
parser.add_argument("--ionizing", type=parsing.int_tuple, help="the range of the FUV luminosity of the ionizing stellar population")
parser.add_argument("--dust", type=parsing.int_tuple, help="the range of the dust mass")

# Logging options
parser.add_argument("--debug", action="store_true", help="add this option to enable debug output")
parser.add_argument("--report", action='store_true', help='write a report file')
parser.add_argument("--config", type=str, help="the name of a configuration file")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Set the modeling path
if arguments.path is None: arguments.path = filesystem.cwd()

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = filesystem.join(arguments.path, time.unique_name("log") + ".txt") if arguments.report else None

# Determine the log level
level = "DEBUG" if arguments.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting explore_advanced ...")

# -----------------------------------------------------------------

# Create an AdvancedParameterExplorer object
explorer = AdvancedParameterExplorer.from_arguments(arguments)

# Run the parameter exploration
explorer.run()

# -----------------------------------------------------------------
