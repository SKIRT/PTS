#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.magic.extract Run galaxy, star and sky extraction on an image
#

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import argparse

# Import astronomical modules
from astropy import log

# Import AstroMagic modules
from astromagic import Image
from astromagic.magic.galaxyextraction import GalaxyExtractor
from astromagic.magic.starextraction import StarExtractor

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument("image", type=str, help="the name of the input image")
parser.add_argument("--regions", action="store_true", help="save regions")
parser.add_argument("--masks", action="store_true", help="save masks")
parser.add_argument("--out", type=str, help="the name of the output directory")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# If an output directory is given
if arguments.out is not None:
    
    # Determine the full path to the output directory
    output_path = os.path.abspath(arguments.out)
    
    # Create the directory if it does not yet exist
    if not os.path.isdir(output_path): os.makedirs(output_path)

# If no output directory is given, place the output in the current working directory
else: output_path = os.getcwd()

# Determine the full path to the image
image_path = os.path.abspath(arguments.image)

# Open the image
image = Image(image_path)

# Create GalaxyExtractor, StarExtractor and SkyExtractor objects
galaxyex = GalaxyExtractor()
starex = StarExtractor()
skyex = SkyExtractor()

# -----------------------------------------------------------------

# Set the appropriate configuration settings for saving the region files
if arguments.regions:
    
    galaxyex.config.save_region = True
    galaxyex.config.saving.region_path = os.path.join(output_path, "galaxies.reg")
    starex.config.save_region = True
    starex.config.saving.region_path = os.path.join(output_path, "stars.reg")

# Set the appropriate configuration settings for saving the masked frames
if arguments.masks:
    
    galaxyex.config.save_masked_frame = True
    galaxyex.config.saving.masked_frame_path = os.path.join(output_path, "masked_galaxies.fits")
    starex.config.save_masked_frame = True
    starex.config.saving.masked_frame_path = os.path.join(output_path, "masked_stars.fits")

# -----------------------------------------------------------------

# Run the galaxy extractor on the primary image frame
galaxyex.run(image.frames.primary)

# Run the star extractor on the primary image frame, passing the galaxy extractor
starex.run(image.frames.primary, galaxyex)

# Run the sky extractor on the primary image frame
skyex.run(image.frames.primary, galaxyex, starex)

# -----------------------------------------------------------------
