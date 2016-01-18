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

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Get the path to the current working directory
working_directory = os.getcwd()

# -----------------------------------------------------------------

# Create a GalaxyModeler object
modeler = GalaxyModeler.from_arguments(arguments)

# Run the modeling
modeler.run(working_directory)

# -----------------------------------------------------------------
