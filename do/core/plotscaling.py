#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.plotscaling Make plots for the scaling benchmark output
#

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import argparse

# Import the relevant PTS classes and modules
from pts.core.plot.scaling import ScalingPlotter

# -----------------------------------------------------------------

# Create the command-line parser and a set of subparsers
parser = argparse.ArgumentParser()
parser.add_argument('--plotfit', action='store_true', help='include the fit to the speedups of a theoretical scaling relation in the plot')

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Create a ScalingPlotter instance
plotter = ScalingPlotter()

# Run the scaling plotter
plotter.run()

# -----------------------------------------------------------------
