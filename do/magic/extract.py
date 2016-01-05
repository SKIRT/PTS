#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.extract Run galaxy, star and sky extraction on an image

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import argparse

# Import the relevant AstroMagic classes and modules
from pts.magic.core import Image
from pts.magic.basics import Mask
from pts.magic import Extractor

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument("image", type=str, help="the name of the input image")
parser.add_argument('--config', type=str, help='the name of a configuration file', default=None)
parser.add_argument("--regions", action="store_true", help="save regions")
parser.add_argument("--masks", action="store_true", help="save masks")
parser.add_argument("-i", "--input", type=str, help="the name of the input directory")
parser.add_argument("-o", "--output", type=str, help="the name of the output directory")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# -- Input --

# If an input directory is given
if arguments.input is not None:

    # Determine the full path to the input directory
    arguments.input_path = os.path.abspath(arguments.input)

    # Give an error if the input directory does not exist
    if not os.path.isdir(arguments.input_path): raise argparse.ArgumentError("The input directory does not exist")

# If no input directory is given, assume the input is placed in the current working directory
else: arguments.input_path = os.getcwd()

# -- Output --

# If an output directory is given
if arguments.output is not None:
    
    # Determine the full path to the output directory
    arguments.output_path = os.path.abspath(arguments.output)
    
    # Create the directory if it does not yet exist
    if not os.path.isdir(arguments.output_path): os.makedirs(arguments.output_path)

# If no output directory is given, place the output in the current working directory
else: arguments.output_path = os.getcwd()

# -----------------------------------------------------------------

# Determine the full path to the image
image_path = os.path.abspath(arguments.image)

# Open the image
image = Image(image_path)

# -----------------------------------------------------------------

# Create a mask for the nans in the primary
mask = Mask.is_nan(image.frames.primary)

# Create an Extractor instance and configure it according to the command-line arguments
extractor = Extractor.from_arguments(arguments)

# Run the extractor
extractor.run(image.frames.primary, mask)

# -----------------------------------------------------------------
