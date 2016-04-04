#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.explore Explore the parameter space for the radiative transfer modeling.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import argparse

# Import the relevant PTS classes and modules
from pts.modeling.fitting.parameterexploration import ParameterExplorer
from pts.core.tools import logging, time, parsing, filesystem

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()

# Basic
parser.add_argument("path", type=str, nargs='?', help="the modeling path")

parser.add_argument("--young_nvalues", type=int, help="the number of different values for the young stellar luminosity")
parser.add_argument("--ionizing_nvalues", type=int, help="the number of different values for the ionizing stellar luminosity")
parser.add_argument("--dust_nvalues", type=int, help="the number of different values for the dust mass")

parser.add_argument("--young_range", type=parsing.float_tuple, help="the range of the luminosity of non-ionizing stars")
parser.add_argument("--ionizing_range", type=parsing.float_tuple, help="the range of the luminosity of ionizing stars")
parser.add_argument("--dust_range", type=parsing.float_tuple, help="the range of the dust mass")

parser.add_argument("--young_log", action="store_true", help="use logarithmic spacing of the young stellar luminosity values")
parser.add_argument("--ionizing_log", action="store_true", help="use logarithmic spacing of the ionizing stellar luminosity values")
parser.add_argument("--dust_log", action="store_true", help="use logarithmic spacing of the dust mass values")

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
logfile_path = filesystem.join(arguments.path, time.unique_name("exploration") + ".txt") if arguments.report else None

# Determine the log level
level = "DEBUG" if arguments.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting parameter exploration ...")

# -----------------------------------------------------------------

# Create a ParameterExplorer object
explorer = ParameterExplorer.from_arguments(arguments)

# Run the parameter exploration
explorer.run()

# -----------------------------------------------------------------
