#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.interpolate Interpolate an image within regions defined by the user.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import argparse

# Import the relevant AstroMagic classes and modules
from pts.magic import ImageImporter
from pts.magic.tools import plotting, interpolation
from pts.magic.core import Frame
from pts.magic.basics import Mask, Region
from pts.core.tools import logging, time, filesystem

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument("image", type=str, help="the name of the input image")
parser.add_argument("region", type=str, help="the name of the region file")
parser.add_argument('--report', action='store_true', help='write a report file')
parser.add_argument("-i", "--input", type=str, help="the name of the input directory")
parser.add_argument("-o", "--output", type=str, help="the name of the output directory")
parser.add_argument("--debug", action="store_true", help="enable debug logging mode")

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

# Determine the log file path
logfile_path = os.path.join(arguments.output_path, time.unique_name("interpolation") + ".txt") if arguments.report else None

# Determine the log level
level = "DEBUG" if arguments.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.info("Starting interpolation script ...")

# -----------------------------------------------------------------

# Inform the user
log.info("Loading the image ...")

# Determine the full path to the image
image_path = os.path.abspath(arguments.image)

# Import the image
importer = ImageImporter()
importer.run(image_path)

# Get the primary image frame
frame = importer.image.frames.primary

# Get the original header
header = importer.image.original_header

# -----------------------------------------------------------------

# Inform the user
log.info("Loading the region ...")

# Load in the region
region_path = filesystem.join(arguments.input_path, arguments.region)
region = Region.from_file(region_path, only=["circle", "ellipse", "polygon"])

# Inform the user
log.info("Creating a mask from the region ...")

# Create a mask from the region
mask = region.to_mask(frame.xsize, frame.ysize)

#plotting.plot_box(mask)

# Inform the user
log.info("Interpolating the frame within the masked pixels ...")

# Create a mask of the pixels that are NaNs
nans = Mask.is_nan(frame)

# Set the NaN pixels to zero in the frame
frame[nans] = 0.0

# Interpolate the frame in the masked pixels
data = interpolation.in_paint(frame, mask)
new_frame = Frame(data)

# Set the original NaN pixels back to NaN
new_frame[nans] = float("nan")

# Inform the user
log.info("Saving the result ...")

# Save the result
path = filesystem.join(arguments.output_path, arguments.image)
new_frame.save(path, header=header)

# -----------------------------------------------------------------
