#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.run Run a SKIRT or FitSKIRT simulation with the best performance, based on the current load of the
#  system.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import argparse

# Import the relevant PTS classes and modules
from pts.core.simulation import SkirtLauncher, FitSkirtLauncher

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument("file", type=str, help="the name of the ski/fski file")
parser.add_argument("--input", type=str, help="the simulation input path")
parser.add_argument("--output", type=str, help="the simulation output path")
parser.add_argument("--remote", action="store_true", help="run the simulation remotely")
parser.add_argument("--relative", action="store_true", help="treats the given input and output paths as being relative to the ski/fski file")
parser.add_argument("--brief", action="store_true", help="enable brief console logging")
parser.add_argument("--verbose", action="store_true", help="enable verbose logging mode")
parser.add_argument("--memory", action="store_true", help="enable memory logging mode")
parser.add_argument("--allocation", action="store_true", help="enable memory (de)allocation logging mode")
parser.add_argument("--emulate", action="store_true", help="emulate the simulation while limiting computation")
parser.add_argument("--extractmemory", action="store_true", help="extract the memory usage from the log files")
parser.add_argument("--plotseds", action="store_true", help="make plots of the output SEDs")
parser.add_argument("--plotgrids", action="store_true", help="make plots of the dust grid")
parser.add_argument("--plotprogress", action="store_true", help="make plots of the progress of the different processes as a function of time")
parser.add_argument("--plottimeline", action="store_true", help="make a plot of the timeline for the different processes")
parser.add_argument("--plotmemory", action="store_true", help="make a plot of the memory consumption as a function of time")
parser.add_argument("--makergb", action="store_true", help="add this option to make RGB images from the SKIRT output")
parser.add_argument("--makewave", action="store_true", help="add this option to make a wave movie from the SKIRT output")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Determine the full path to the parameter file
arguments.filepath = os.path.abspath(arguments.file)

# Determine the full path to the input and output directories
if arguments.input is not None: arguments.input_path = os.path.abspath(arguments.input)
if arguments.output is not None: arguments.output_path = os.path.abspath(arguments.output)

# -----------------------------------------------------------------

# If the parameter file describes a SKIRT simulation
if arguments.filepath.endswith(".ski"):

    # Create a SkirtLauncher instance and run it
    launcher = SkirtLauncher.from_arguments(arguments)
    launcher.run()

# If the parameter file describes a FitSKIRT simulation
elif arguments.filepath.endswith(".fski"):

    # Create a FitSkirtLauncher instance and run it
    launcher = FitSkirtLauncher.from_arguments(arguments)
    launcher.run()

# If the parameter file has a different extension
else: raise argparse.ArgumentError("The parameter file is not a ski or fski file")

# -----------------------------------------------------------------
