#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.preparedata Do the preparation step of the galaxy modeling

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import argparse

# Import the relevant PTS classes and modules
from pts.modeling.core import DataPreparer

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument("image", type=str, help="the name of the image for which to run the preparation")
parser.add_argument("--reference", type=str, help="the name of the reference image")
parser.add_argument("--config", type=str, help="the name of a configuration file", default=None)

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Create a DataPreparer instance
preparer = DataPreparer.from_arguments(arguments)

# Run the data preparation
preparer.run()

# -----------------------------------------------------------------
