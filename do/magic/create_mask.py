#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.create_mask Create a mask from a region file, for a particular image

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import argparse

# Import the relevant PTS classes and modules
from pts.magic.core.frame import Frame
from pts.magic.region.list import load_as_pixel_region_list
from pts.magic.basics.mask import Mask

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument("region", type=str, help="the name of the region file")
parser.add_argument("image", type=str, help="the name of the image file")
parser.add_argument("value", type=float, nargs='?', help="the fill value", default='nan')
parser.add_argument("--data", action="store_true", help="use the original data for pixels that are not masked")
parser.add_argument('--invert', action="store_true", help="invert the mask so that mask covers outside regions")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Load the image
frame = Frame.from_file(arguments.image)

# Load the region
region_name = os.path.splitext(os.path.basename(arguments.region))[0]
region = load_as_pixel_region_list(arguments.region, frame.wcs)

# Create the mask
mask = region.to_mask(x_size=frame.xsize, y_size=frame.ysize)

# Calculate the inverse, if requested
if arguments.invert: mask = mask.inverse()

# -----------------------------------------------------------------

if arguments.data:
    
    new_frame = frame
    frame[mask] = arguments.value
    
# Create a frame for the total mask
else: new_frame = Frame(mask.astype(int))

# Write out the new frame
new_frame.saveto(region_name + ".fits")

# -----------------------------------------------------------------
