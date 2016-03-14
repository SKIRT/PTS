#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.subtract Run sky subtraction on an astronmical image.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import argparse

# Import the relevant PTS classes and modules
from pts.magic.misc.imageimporter import ImageImporter
from pts.magic.sky.skysubtractor import SkySubtractor
from pts.magic.basics.mask import Mask
from pts.core.tools import configuration

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()

# Basic
parser.add_argument("image", type=str, help="the name of the input image")

# Logging
parser.add_argument("--debug", action="store_true", help="enable debug logging mode")
parser.add_argument('--report', action='store_true', help='write a report file')


parser.add_argument("mask", type=str, help="the name of the mask image resulting from the extraction procedure")
parser.add_argument("galaxies", type=str, help="the name of the file specifying the galaxy regions")

#parser.add_argument("saturation", type=str, nargs='?', help="the name of the file specifying the saturation regions", default=None)

parser.add_argument('--config', type=str, help='the name of a configuration file')
parser.add_argument("--settings", type=configuration.from_string, help="settings")
parser.add_argument("-i", "--input", type=str, help="the name of the input directory")
parser.add_argument("-o", "--output", type=str, help="the name of the output directory")

parser.add_argument("--bad", type=str, help="the name of the file specifying regions that have to be added to the mask of bad pixels")


# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# -- Input --

# If an input directory is given
if arguments.input is not None:

    # Determine the full path to the input directory
    arguments.input_path = os.path.abspath(arguments.input)

    # Give an error if the input directory does not exist
    if not os.path.isdir(arguments.input_path): raise argparse.ArgumentError(arguments.input_path, "The input directory does not exist")

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

# Determine the full path to the bad region file
bad_region_path = os.path.join(arguments.input_path, arguments.bad) if arguments.bad is not None else None

# Import the image
importer = ImageImporter()
importer.run(image_path, bad_region_path=bad_region_path)

# Determine the full path to the mask
mask_path = os.path.abspath(arguments.mask)

# Open the mask frame
mask = Mask.from_file(mask_path)

# -----------------------------------------------------------------

# Determine the full path to the galaxy region file and the saturation region file
galaxy_region_path = os.path.join(arguments.input_path, arguments.galaxies)
saturation_region_path = os.path.join(arguments.input_path, arguments.saturation) if arguments.saturation is not None else None

# Create a SkySubtractor instance and configure it according to the command-line arguments
subtractor = SkySubtractor.from_arguments(arguments)

# Run the subtractor
subtractor.run(importer.image.frames.primary, mask, galaxy_region_path, saturation_region_path, bad_mask=importer.mask)

# Save the result
subtractor.write_result(importer.image.original_header)

# -----------------------------------------------------------------
