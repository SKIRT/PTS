#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.magic.subtract Run sky subtraction on an image
#

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import argparse

# Import the relevant AstroMagic classes and modules
from pts.magic.core import Image, Frame
from pts.magic import SkySubtractor

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument("image", type=str, help="the name of the input image")
parser.add_argument("mask", type=str, help="the name of the mask image")
parser.add_argument('--config', type=str, help='the name of a configuration file', default=None)
parser.add_argument("--out", type=str, help="the name of the output directory")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# If an output directory is given
if arguments.out is not None:
    
    # Determine the full path to the output directory
    arguments.output_path = os.path.abspath(arguments.out)
    
    # Create the directory if it does not yet exist
    if not os.path.isdir(arguments.output_path): os.makedirs(arguments.output_path)

# If no output directory is given, place the output in the current working directory
else: arguments.output_path = os.getcwd()

# -----------------------------------------------------------------

# Determine the full path to the image
image_path = os.path.abspath(arguments.image)

# Open the image
image = Image(image_path)

# Determine the full path to the mask
mask_path = os.path.abspath(arguments.mask)

# Open the mask frame
mask = Frame.from_file(mask_path)

# -----------------------------------------------------------------

# Create a SkySubtractor instance and configure it according to the command-line arguments
subtractor = SkySubtractor.from_arguments(arguments)

# Run the subtractor
subtractor.run(image.frames.primary, mask)

# -----------------------------------------------------------------
