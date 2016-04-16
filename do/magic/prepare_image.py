#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.prepare_image Prepare an image with PTS.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import argparse

# Import the relevant PTS classes and modules
from pts.magic.prepare.imagepreparation import ImagePreparer
from pts.core.tools import logging, time, filesystem
from pts.magic.core.image import Image
from pts.magic.basics.region import Region

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()

# Basic
parser.add_argument("image", type=str, nargs='?', help="the name/path of the image for which to run the preparation")
parser.add_argument("reference", type=str, help="the name/path of the reference image (to which the image is rebinned)")
parser.add_argument("kernel", type=str, help="the name/path of the kernel file for the convolution")

# Input and output
parser.add_argument("--input", type=str, help="the input path (output of find_sources step)")
parser.add_argument("--output", type=str, help="the output path")

# Logging
parser.add_argument("--debug", action="store_true", help="enable debug logging mode")
parser.add_argument("--report", action='store_true', help="write a report file")

parser.add_argument("--steps", action="store_true", help="write the results of intermediate steps")
parser.add_argument("--config", type=str, help="the name of a configuration file")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Determine the full input and output paths
if arguments.output is None: arguments.output = filesystem.cwd()
if arguments.input is None: arguments.input = filesystem.cwd()
arguments.input = filesystem.absolute(arguments.input)
arguments.output = filesystem.absolute(arguments.output)

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = filesystem.join(arguments.output, time.unique_name("log") + ".txt") if arguments.report else None

# Determine the log level
level = "DEBUG" if arguments.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting prepare_image ...")

# -----------------------------------------------------------------

# Determine the path to the input image
image_path = filesystem.absolute(arguments.image)

# Load the image
image = Image.from_file(image_path)

# Determine the absolute path to the reference image
arguments.reference = filesystem.absolute(arguments.reference)

# Determine the absolute path to the convolution kernel
arguments.kernel = filesystem.absolute(arguments.kernel)

# Determine the path to the galaxy region
galaxy_region_path = filesystem.join(arguments.input, "galaxies.reg")

# Load the galaxy region
galaxy_region = Region.from_file(galaxy_region_path)

# Determine the path to the star region
star_region_path = filesystem.join(arguments.input, "stars.reg")

# Load the star region
star_region = Region.from_file(star_region_path)

# Determine the path to the saturation region
saturation_region_path = filesystem.join(arguments.input, "saturation.reg")

# Load the saturation region
saturation_region = Region.from_file(saturation_region_path)

# Determine the path to the region of other sources
other_region_path = filesystem.join(arguments.input, "other_sources.reg")

# Load the region of other sources
other_region = Region.from_file(other_region_path)

# Load the image with segmentation maps
segments_path = filesystem.join(arguments.input, "segments.fits")
segments = Image.from_file(segments_path)

# Get the segmentation maps
galaxy_segments = segments.frames.galaxies
star_segments = segments.frames.stars
other_segments = segments.frames.other_sources

# -----------------------------------------------------------------

# Create an ImagePreparer instance
preparer = ImagePreparer.from_arguments(arguments)

# Run the image preparation
preparer.run(image, galaxy_region, star_region, saturation_region, other_region, galaxy_segments, star_segments, other_segments)

# -----------------------------------------------------------------
