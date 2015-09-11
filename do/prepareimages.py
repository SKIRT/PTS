#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.prepareimages Prepare images for SKIRT radiative transfer simulations
#

# *****************************************************************

# Import standard modules
import os.path
import argparse

# Import relevant PTS modules
from pts.imagepreparation import ImagePreparation

# *****************************************************************

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument('filter', type=str, help='the filter for which to run the data preparation', nargs='?', default=None)
parser.add_argument("--stage", type=str, help="the preparation stage")
parser.add_argument('--plot', action='store_true', help='plot the results of intermediate steps')
parser.add_argument('--save', action='store_true', help='save intermediate results')

# Parse the command line arguments
args = parser.parse_args()

# Set the command-line options
filter = args.filter
stage = args.stage
plot = args.plot
save = args.save

# *****************************************************************

# Get the path to the current working directory
working_directory = os.getcwd()

# *****************************************************************

# Create a ImagePreparation object
preparation = ImagePreparation(working_directory, filter, plot, save)

# Run the image preparation
if stage is None: preparation.run()
elif stage == "prepare": preparation.prepare()
elif stage == "galfit": preparation.bulge_and_disk()
elif stage == "maps": preparation.make_maps()
else: raise ValueError("Unkown stage")

# *****************************************************************
