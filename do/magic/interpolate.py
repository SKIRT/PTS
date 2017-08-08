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

# Import the relevant PTS classes and modules
from pts.magic.misc.imageimporter import ImageImporter
from pts.magic.tools import interpolation
from pts.magic.core.frame import Frame
from pts.magic.basics.mask import Mask
from pts.magic.region.list import PixelRegionList
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.basics.log import log

# -----------------------------------------------------------------

default_method = "biharmonic"

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()
definition.add_required("image", "string", "name of the input image")
definition.add_required("regions", "string", "name of the regions file over which to interpolate")
definition.add_flag("mask", "write out the mask")
definition.add_optional("color", "string", "only interpolate over the shapes with this color")
definition.add_optional("ignore_color", "string", "ignore shapes with this particular color")
definition.add_optional("shapes", "string_list", "only interpolate over these kinds of shapes")
definition.add_optional("input", "string", "name of the input directory", letter="i")
definition.add_optional("output", "string", "name of the output directory", letter="o")
definition.add_optional("method", "string", "interpolation method to use", default=default_method)
config = parse_arguments("interpolate", definition)

# -----------------------------------------------------------------

# If an input directory is given
if config.input is not None:

    # Determine the full path to the input directory
    input_path = fs.absolute_path(config.input)

    # Give an error if the input directory does not exist
    if not fs.is_directory(input_path): raise IOError(input_path, "The input directory does not exist")

# If no input directory is given, assume the input is placed in the current working directory
else: input_path = fs.cwd()

# -----------------------------------------------------------------

# If an output directory is given
if config.output is not None:

    # Determine the full path to the output directory
    output_path = fs.absolute_path(config.output)
    
    # Create the directory if it does not yet exist
    if not fs.is_directory(output_path): fs.create_directory(output_path)

# If no output directory is given, place the output in the current working directory
else: output_path = fs.cwd()

# -----------------------------------------------------------------

# Inform the user
log.info("Loading the image ...")

# Determine the full path to the image
image_path = fs.absolute_path(config.image)

# Import the image
importer = ImageImporter()
importer.run(image_path)

# Get the primary image frame
frame = importer.image.primary

# Get the original header
header = importer.image.original_header

# -----------------------------------------------------------------

# Inform the user
log.info("Loading the region ...")

# Load in the region
region_path = fs.join(input_path, config.region)
region = PixelRegionList.from_file(region_path, only=config.shapes, color=config.color, ignore_color=ignore.ignore_color)

# Inform the user
log.info("Creating a mask from the region ...")

# Create a mask from the region
mask = region.to_mask(frame.xsize, frame.ysize)

# Inform the user
log.info("Interpolating the frame within the masked pixels ...")

# Create a mask of the pixels that are NaNs
nans = Mask.is_nan(frame)

# Set the NaN pixels to zero in the frame
frame[nans] = 0.0

# Interpolate the frame in the masked pixels
if config.method == "biharmonic": data = interpolation.inpaint_biharmonic(frame, mask)
elif config.method == "local_mean": data = interpolation.in_paint(frame, mask, method="localmean")
else: raise ValueError("Invalid interpolation method (should be 'biharmonic' or 'local_mean')")
new_frame = Frame(data)

# Set the original NaN pixels back to NaN
new_frame[nans] = float("nan")

# Inform the user
log.info("Saving the result ...")

# Save the result
path = fs.join(output_path, config.image)
new_frame.saveto(path, header=header)

# Write the mask
if config.mask:

    path = fs.join(output_path, "mask.fits")
    new_frame[mask] = float('nan')
    new_frame.saveto(path, header=header)

# -----------------------------------------------------------------
