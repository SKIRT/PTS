#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.modelgalaxy Model a galaxy with Astromagic and SKIRT

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import argparse

# Import the relevant PTS classes and modules
from pts.modeling import GalaxyModeler

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument("--config", type=str, help='the name of a configuration file', default=None)
parser.add_argument("--stage", type=str, help='the preparation stage')
parser.add_argument('--report', action='store_true', help='write a report file')
parser.add_argument('--plot', action='store_true', help='plot the result of intermediate steps')
parser.add_argument("--debug", action="store_true", help="enable debug logging mode")
parser.add_argument("--steps", action="store_true", help="write out the results of intermediate steps")
parser.add_argument("-i", "--input", type=str, help="the name of the input directory")
parser.add_argument("-o", "--output", type=str, help="the name of the output directory")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# -- Input --

# If an input directory is given
if arguments.input is not None:

    # Determine the full path to the input directory
    arguments.input_path = os.path.abspath(arguments.input)

    # Give an error if the input directory does not exist
    if not os.path.isdir(arguments.input_path): raise argparse.ArgumentError(arguments.input_path, "The input directory does not exist")

# If no input directory is given, assume the input is placed in the current working directory
else: arguments.input_path = os.getcwd()

# -- Output --

# If an output directory is given
if arguments.output is not None:

    # Determine the full path to the output directory
    arguments.output_path = os.path.abspath(arguments.output)

    # Create the directory if it does not yet exist
    if not os.path.isdir(arguments.output_path): os.makedirs(arguments.output_path)

# If no output directory is given, place the output in the current working directory
else: arguments.output_path = os.getcwd()

# Get the path to the current working directory
working_directory = os.getcwd()

# -----------------------------------------------------------------

# Create a GalaxyModeler object
modeler = GalaxyModeler.from_arguments(arguments)

# Run the modeling
modeler.run(working_directory)

# -----------------------------------------------------------------
