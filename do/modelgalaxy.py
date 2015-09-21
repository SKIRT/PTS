#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.modelgalaxy Model a galaxy with Astromagic and SKIRT
#

# *****************************************************************

# Import standard modules
import os.path
import argparse

# Import relevant PTS modules
from pts.galaxymodeler import GalaxyModeler

# *****************************************************************

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument('--image', type=str, help='provide this argument to only prepare one specific image')
parser.add_argument('--stage', type=str, help='the preparation stage')
parser.add_argument('--plot', action='store_true', help='plot the results of intermediate steps')
parser.add_argument('--save', action='store_true', help='save intermediate results')

# Parse the command line arguments
args = parser.parse_args()

# Set the command-line options
filter_name = args.filter
stage = args.stage
plot = args.plot
save = args.save

# *****************************************************************

# Get the path to the current working directory
working_directory = os.getcwd()

# *****************************************************************

# Create a GalaxyModeler object
modeler = GalaxyModeler(working_directory, filter_name, plot, save)

# Run the modeling procedure
if stage is None: modeler.run()
elif stage == "prepare": modeler.prepare_images()
elif stage == "galfit": modeler.fit_bulge_and_disk()
elif stage == "maps": modeler.make_maps()
elif stage == "fit": modeler.fit_sed()
else: raise ValueError("Unkown stage")

# *****************************************************************
