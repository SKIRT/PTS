#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.photometry Calculate the photometry of the different input images

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import argparse

# Import the relevant PTS classes and modules
from pts.modeling.core import PhotoMeter

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument("image", type=str, help="the name of the image for which to calculate the photometry")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Create a PhotoMeter object
photometer = PhotoMeter.from_arguments(arguments)

# Run the photometry
photometer.run()

# -----------------------------------------------------------------
