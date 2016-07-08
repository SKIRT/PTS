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
from pts.magic.prepare.preparer import ImagePreparer
from pts.core.tools import logging, time, tables, parsing
from pts.core.tools import filesystem as fs
from pts.magic.core.image import Image
from pts.magic.basics.region import Region
from pts.magic.misc.calibration import CalibrationError
from pts.magic.misc.extinction import GalacticExtinction
from pts.core.basics.filter import Filter
from pts.magic.misc.kernels import AnianoKernels, aniano_names, variable_fwhms
from pts.core.basics.configuration import Configuration

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()

# Basic
parser.add_argument("image", type=str, help="the name/path of the image for which to run the preparation")
parser.add_argument("convolve_to", type=str, help="the name of the band to convolve the image to")
parser.add_argument("rebin_to", type=str, nargs='?', help="the name/path of the reference image to which the image is rebinned")

# Advanced options
parser.add_argument("--sky_annulus_outer", type=float, help="the factor to which the ellipse describing the principal galaxy should be multiplied to represent the outer edge of the sky annulus")
parser.add_argument("--sky_annulus_inner", type=float, help="the factor to which the ellipse describing the principal galaxy should be multiplied to represent the inner edge of the sky annulus")
parser.add_argument("--convolution_remote", type=str, help="the name of the remote host to be used for the convolution step")
parser.add_argument("--sky_region", type=str, help="the name/path of a file with manually selected regions for the sky estimation (not apertures but extended regions of any shape and number) (in sky coordinates!)")
parser.add_argument("--error_frames", type=parsing.string_list, help="the names of planes in the input image which have to be regarded as error maps (seperated by commas)")

# Input and output
parser.add_argument("--input", type=str, help="the input path (output of find_sources step)")
parser.add_argument("--output", type=str, help="the output path")

# Logging
parser.add_argument("--debug", action="store_true", help="enable debug logging mode")
parser.add_argument("--report", action='store_true', help="write a report file")

parser.add_argument("--steps", action="store_true", help="write the results of intermediate steps")
parser.add_argument("--config", type=str, help="the name of a configuration file")

# Visualisation
parser.add_argument("--visualise", action="store_true", help="make visualisations")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Determine the full input and output paths
if arguments.output is None: arguments.output = fs.cwd()
if arguments.input is None: arguments.input = fs.cwd()
arguments.input = fs.absolute(arguments.input)
arguments.output = fs.absolute(arguments.output)

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = fs.join(arguments.output, time.unique_name("log") + ".txt") if arguments.report else None

# Determine the log level
level = "DEBUG" if arguments.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting prepare_image ...")

# -----------------------------------------------------------------

# Determine the path to the input image
image_path = fs.absolute(arguments.image)

# Load the image
image = Image.from_file(image_path)

# -----------------------------------------------------------------

# Inform the user
log.info("Loading regions ...")

# Determine the path to the galaxy region
galaxy_region_path = fs.join(arguments.input, "galaxies.reg")

# Load the galaxy region
galaxy_region = Region.from_file(galaxy_region_path)

# Determine the path to the star region
star_region_path = fs.join(arguments.input, "stars.reg")

# Load the star region
star_region = Region.from_file(star_region_path) if fs.is_file(star_region_path) else None

# Determine the path to the saturation region
saturation_region_path = fs.join(arguments.input, "saturation.reg")

# Load the saturation region
saturation_region = Region.from_file(saturation_region_path) if fs.is_file(saturation_region_path) else None

# Determine the path to the region of other sources
other_region_path = fs.join(arguments.input, "other_sources.reg")

# Load the region of other sources
other_region = Region.from_file(other_region_path) if fs.is_file(other_region_path) else None

# Inform the user
log.debug("Loading segmentation frames ...")

# Load the image with segmentation maps
segments_path = fs.join(arguments.input, "segments.fits")
segments = Image.from_file(segments_path, no_filter=True)

# Get the segmentation maps
galaxy_segments = segments.frames.galaxies
star_segments = segments.frames.stars
other_segments = segments.frames.other_sources

# Load the statistics file
statistics_path = fs.join(arguments.input, "statistics.dat")

# Inform the user
log.debug("Loading the FWHM ...")

# Get the FWHM from the statistics file
fwhm = None
with open(statistics_path) as statistics_file:
    for line in statistics_file:
        if "FWHM" in line: fwhm = parsing.get_quantity(line.split("FWHM: ")[1].replace("\n", ""))

# -----------------------------------------------------------------

# Get the center coordinate of the frame
center_coordinate = image.coordinate_range[0]

# Get the filter name
if image.filter is None: raise RuntimeError("Filter not recognized!")
filter_name = str(image.filter)

# -----------------------------------------------------------------

# Debugging
log.debug("Getting galactic extinction ...")

# Get the galactic extinction for this image
arguments.attenuation = GalacticExtinction(center_coordinate).extinction_for_filter(image.filter)

# -----------------------------------------------------------------

# Get the calibration error
arguments.calibration = CalibrationError.from_filter(image.filter)

# -----------------------------------------------------------------

# If visualisation is enabled, set the visualisation path (=output path)
if arguments.visualise: visualisation_path = arguments.output
else: visualisation_path = None

# -----------------------------------------------------------------

# Inform the user
log.info("Looking up the necessary kernel file ...")

# Get the filter to which to convolve to
convolve_to_filter = Filter.from_string(arguments.convolve_to)

# Create an AnianoKernels instance
kernels = AnianoKernels()

# Get the path to the appropriate convolution kernel
kernel_path = kernels.get_kernel_path(image.filter, convolve_to_filter, fwhm=fwhm)

# Set the kernel path
arguments.kernel = kernel_path

# -----------------------------------------------------------------

# Determine the absolute path to the reference image
arguments.rebin_to = fs.absolute(arguments.rebin_to)

# Determine the full path to the sky region file
if arguments.sky_region is not None: arguments.sky_region = fs.absolute(arguments.sky_region)

# -----------------------------------------------------------------

# Create an ImagePreparer instance
preparer = ImagePreparer.from_arguments(arguments)

# Run the image preparation
preparer.run(image, galaxy_region, star_region, saturation_region, other_region, galaxy_segments, star_segments, other_segments, visualisation_path)

# -----------------------------------------------------------------
