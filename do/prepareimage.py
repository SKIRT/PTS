#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.prepareimage Prepare an image with Astromagic
#

# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import argparse

# Import the relevant PTS modules
from modeling.imagepreparation import ImagePreparation

# Import Astromagic modules
from astromagic import Image

# *****************************************************************

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument('path', type=str, help='the path to the image')
parser.add_argument('--config', type=str, help='the name of a configuration file', default=None)

# Parse the command line arguments
args = parser.parse_args()

# Set the command-line options
path = args.path
config_file = args.config

# *****************************************************************

# Create an image
image = Image(path)

# Create a ImagePreparation object
prep = ImagePreparation(config_file)

# Run the image preparation on the image
prep.run(image)

# *****************************************************************
