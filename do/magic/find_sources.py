#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.find_sources Find sources (galaxies, stars) in an astronomical image.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import argparse

# Import the relevant PTS classes and modules
from pts.magic.misc.imageimporter import ImageImporter
from pts.magic.sources.finder import SourceFinder
from pts.magic.catalog.importer import CatalogImporter
from pts.magic.core.image import Image
from pts.magic.basics.region import Region
from pts.core.tools import configuration
from pts.core.tools import logging, time, filesystem

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument("image", type=str, help="the name of the input image")
parser.add_argument("--debug", action="store_true", help="enable debug logging mode")
parser.add_argument('--report', action='store_true', help='write a report file')

parser.add_argument('--config', type=str, help='the name of a configuration file', default=None)
parser.add_argument("--settings", type=configuration.from_string, help="settings")

parser.add_argument("-i", "--input", type=str, help="the name of the input directory")
parser.add_argument("-o", "--output", type=str, help="the name of the output directory")

parser.add_argument("--synchronize", action="store_true", help="synchronize with DustPedia catalog")
parser.add_argument("--filecatalog", action="store_true", help="use file catalogs")
parser.add_argument("--interpolation_method", type=str, help="the interpolation method to use")

parser.add_argument("--ignore", type=str, help="the name of the file specifying regions to ignore")
parser.add_argument("--special", type=str, help="the name of the file specifying regions with objects needing special attention")
parser.add_argument("--bad", type=str, help="the name of the file specifying regions that have to be added to the mask of bad pixels")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# If an input directory is given
if arguments.input is not None:

    # Determine the full path to the input directory
    input_path = filesystem.absolute(arguments.input)

    # Give an error if the input directory does not exist
    if not filesystem.is_directory(input_path):
        raise argparse.ArgumentError(input_path, "The input directory does not exist")

# If no input directory is given, assume the input is placed in the current working directory
else: input_path = filesystem.cwd()

# -----------------------------------------------------------------

# If an output directory is given
if arguments.output is not None:
    
    # Determine the full path to the output directory
    output_path = filesystem.absolute(arguments.output)
    
    # Create the directory if it does not yet exist
    if not filesystem.is_directory(output_path): filesystem.create_directory(output_path)

# If no output directory is given, place the output in the current working directory
else: output_path = filesystem.cwd()

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = filesystem.join(output_path, time.unique_name("sourcefinder") + ".txt") if arguments.report else None

# Determine the log level
level = "DEBUG" if arguments.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
logging.log.info("Starting find_sources ...")

# -----------------------------------------------------------------

# Determine the full path to the image
image_path = filesystem.absolute(arguments.image)

# Determine the full path to the bad region file
bad_region_path = filesystem.join(input_path, arguments.bad) if arguments.bad is not None else None

# Import the image
importer = ImageImporter()
importer.run(image_path, bad_region_path=bad_region_path)

# -----------------------------------------------------------------

# Create a CatalogImporter instance
catalog_importer = CatalogImporter()

# Run the catalog importer
catalog_importer.run(importer.image.frames.primary) # work with coordinate box instead ? image.coordinate_box ?

# -----------------------------------------------------------------

# If a special region is defined
if arguments.special is not None:

    # Determine the full path to the special region file
    path = filesystem.join(input_path, arguments.special)

    # Inform the user
    log.info("Creating mask covering objects that require special attention from " + path + " ...")

    # Load the region and create a mask from it
    special_region = Region.from_file(path)

# No special region
else: special_region = None

# -----------------------------------------------------------------

# If an ignore region is defined
if arguments.ignore is not None:

    # Determine the full path to the ignore region file
    path = filesystem.join(input_path, arguments.ignore)

    # Inform the user
    log.info("Creating mask covering objects that should be ignored from " + path + " ...")

    # Load the region and create a mask from it
    ignore_region = Region.from_file(path)

# No ignore region
else: ignore_region = None

# -----------------------------------------------------------------

# Create a SourceFinder instance
finder = SourceFinder.from_arguments(arguments)

# Run the source finder
finder.run(importer.image, catalog_importer.galactic_catalog, catalog_importer.stellar_catalog, special_region, ignore_region)

# -----------------------------------------------------------------

# Save the galaxy region
galaxy_region = finder.galaxy_region
path = filesystem.join(output_path, "galaxies.reg")
galaxy_region.save(path)

# Save the star region
star_region = finder.star_region
path = filesystem.join(output_path, "stars.reg")
star_region.save(path)

# Save the saturation region
saturation_region = finder.saturation_region
path = filesystem.join(output_path, "saturation.reg")
saturation_region.save(path)

# -----------------------------------------------------------------

# Create an image with the segmentation maps
segments = Image("segments")

# Add the segmentation map of the galaxies
segments.add_frame(finder.galaxy_segments, "galaxies")

# Add the segmentation map of the saturated stars
segments.add_frame(finder.star_segments, "saturation")

# Add the segmentation map of the other sources
segments.add_frame(finder.other_segments, "other_sources")

# Save the FITS file with the segmentation maps
path = filesystem.join(output_path, "segments.fits")
segments.save(path)

# -----------------------------------------------------------------
