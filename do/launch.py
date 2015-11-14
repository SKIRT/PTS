#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.run Run a SKIRT or FitSKIRT simulation with the best performance, based on the current load of the
#  system.

# -----------------------------------------------------------------

# Import standard modules
import os
import argparse

# Import the relevant PTS modules
from pts.launcher import SkirtLauncher, FitSkirtLauncher, SkirtMemoryLauncher

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument("file", type=str, help="the name of the ski/fski file")
parser.add_argument("inpath", type=str, help="the simulation input path")
parser.add_argument("outpath", type=str, help="the simulation output path")
parser.add_argument("--memory", action="store_true", help="run the memory console application")
parser.add_argument("--relative", action="treats the given input and output paths as being relative to the ski/fski file")
parser.add_argument('--brief', action='store_true', help="enable brief console logging")
parser.add_argument('--verbose', action='store_true', help="enable verbose logging mode for SKIRT")
parser.add_argument('--plotseds', action='store_true', help="make plots of the output SEDs")
parser.add_argument('--plotgrids', action='store_true', help="make plots of the dust grid")
parser.add_argument('--plotprogress', action='store_true', help="make plots of the progress of the different processes as a function of time")
parser.add_argument('--plottimeline', action='store_true', help="make a plot of the timeline for the different processes")
parser.add_argument('--makergb', action='store_true', help="add this option to make RGB images from the SKIRT output")
parser.add_argument('--makewave', action='store_true', help="add this option to make a wave movie from the SKIRT output")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Determine the full path to the parameter file
arguments.filepath = os.path.abspath(arguments.file)

# Determine the full path to the input and output directories
arguments.input_path = os.path.abspath(arguments.inpath)
arguments.output_path = os.path.abspath(arguments.outpath)

exit()

# If the parameter file describes a SKIRT simulation
if arguments.filepath.endswith(".ski"):

    # Check whether the 'memory' flag was provided
    if arguments.memory:

        # Create a SkirtMemoryLauncher instance and run it
        launcher = SkirtMemoryLauncher.from_arguments(arguments)
        launcher.run(filepath, input_path, output_path)

    else:

        # Create a SkirtLauncher instance and run it
        launcher = SkirtLauncher.from_arguments(arguments)
        launcher.run()

# If the parameter file describes a FitSKIRT simulation
elif arguments.filepath.endswith(".fski"):

    # Check whether the 'memory' flag was provided
    if arguments.memory: raise ValueError("The memory console application for FitSKIRT does not exist yet")

    else:

        # Create a FitSkirtLauncher instance and run it
        launcher = FitSkirtLauncher.from_arguments(arguments)
        launcher.run()

# If the parameter file has a different extension
else: raise argparse.ArgumentError("The parameter file is not a ski or fski file")

# -----------------------------------------------------------------
