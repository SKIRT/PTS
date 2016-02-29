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
from pts.modeling.fitting import ParameterExplorer

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Create a ParameterExplorer object
explorer = ParameterExplorer.from_arguments(arguments)

# Run the parameter exploration
explorer.run()

# -----------------------------------------------------------------
