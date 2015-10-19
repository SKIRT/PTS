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
from modeling.galaxymodeler import GalaxyModeler

# *****************************************************************

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument('--image', type=str, help='provide this argument to only prepare one specific image')
parser.add_argument('--stage', type=str, help='the preparation stage')
parser.add_argument('--config', type=str, help='the name of a configuration file', default=None)
parser.add_argument('--plot', action='store_true', help='plot the results of intermediate steps')
parser.add_argument('--save', action='store_true', help='save intermediate results')

# Parse the command line arguments
args = parser.parse_args()

# Set the command-line options
filter_name = args.image
stage = args.stage
plot = args.plot
save = args.save
config_file = args.config

# *****************************************************************

# Get the path to the current working directory
working_directory = os.getcwd()

# *****************************************************************

# Create a GalaxyModeler object
modeler = GalaxyModeler(working_directory, config_file)

# Run the modeling procedure
if stage is None: modeler.run()
elif stage == "preparation":

    # Set the name of the image that needs to be prepared and run the image preparation
    modeler.config.preparation.filter_name = filter_name
    modeler.prepare_images()

elif stage == "decomposition": modeler.decompose()
elif stage == "mapmaking": modeler.make_maps()
elif stage == "fitting": modeler.fit_sed()
else: raise ValueError("Unkown stage")

# *****************************************************************
