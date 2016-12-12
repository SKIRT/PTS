#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.find_and_extract Find and extract sources from an astronomical image.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import argparse

# Import the relevant PTS classes and modules
from pts.magic.core.image import Image
from pts.magic.misc.imageimporter import ImageImporter
from pts.magic.sources.finder import SourceFinder
from pts.magic.catalog.importer import CatalogImporter
from pts.magic.sources.extractor import SourceExtractor
from pts.core.tools import configuration
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs
from pts.magic.region.list import PixelRegionList
from pts.magic.view import MagicViewer

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()

# Configuration
parser.add_argument('--config', type=str, help='the name of a configuration file', default=None)
parser.add_argument("--settings", type=configuration.from_string, help="settings")

# Basic
parser.add_argument("image", type=str, help="the name of the input image")

# Logging options
parser.add_argument("--debug", action="store_true", help="enable debug logging mode")
parser.add_argument("--report", action='store_true', help='write a report file')

parser.add_argument("--synchronize", action="store_true", help="synchronize with DustPedia catalog")

parser.add_argument("--filecatalog", action="store_true", help="use file catalogs")

# Input and output
parser.add_argument("-i", "--input", type=str, help="the name of the input directory")
parser.add_argument("-o", "--output", type=str, help="the name of the output directory")

# Regions
parser.add_argument("--special", type=str, help="the name of the file specifying regions with objects needing special attention")
parser.add_argument("--bad", type=str, help="the name of the file specifying regions that have to be added to the mask of bad pixels")

# Interactive
parser.add_argument("--interactive", action="store_true", help="enable interactive mode")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# If an input directory is given
if arguments.input is not None:

    # Determine the full path to the input directory
    input_path = fs.absolute_path(arguments.input)

    # Give an error if the input directory does not exist
    if not fs.is_directory(input_path): raise argparse.ArgumentError(input_path, "The input directory does not exist")

# If no input directory is given, assume the input is placed in the current working directory
else: input_path = fs.cwd()

# -----------------------------------------------------------------

# If an output directory is given
if arguments.output is not None:
    
    # Determine the full path to the output directory
    output_path = fs.absolute_path(arguments.output)
    
    # Create the directory if it does not yet exist
    if not fs.is_directory(output_path): fs.create_directory(output_path)

# If no output directory is given, place the output in the current working directory
else: output_path = fs.cwd()

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = fs.join(output_path, time.unique_name("log") + ".txt") if arguments.report else None

# Determine the log level
level = "DEBUG" if arguments.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting find_and_extract ...")

# -----------------------------------------------------------------

# Determine the full path to the image
image_path = fs.absolute_path(arguments.image)

# Determine the full path to the bad region file
bad_region_path = fs.join(input_path, arguments.bad) if arguments.bad is not None else None

# Import the image
importer = ImageImporter()
importer.run(image_path, bad_region_path=bad_region_path)

# Get the imported image
image = importer.image

# Get the mask of bad pixels
bad_mask = image.masks.bad if "bad" in image.masks else None

# -----------------------------------------------------------------

# Inform the user
log.info("Importing the galactic and stellar catalogs ...")

# Create a CatalogImporter instance
catalog_importer = CatalogImporter()

# Run the catalog importer
catalog_importer.run(image.frames.primary)

# -----------------------------------------------------------------

# If a special region is defined
if arguments.special is not None:

    # Determine the full path to the special region file
    path = fs.join(input_path, arguments.special)

    # Inform the user
    log.info("Creating mask covering objects that require special attention from " + path + " ...")

    # Load the region and create a mask from it
    special_region = PixelRegionList.from_file(path)

# No special region
else: special_region = None

# -----------------------------------------------------------------

# If an ignore region is defined
if arguments.ignore is not None:

    # Determine the full path to the ignore region file
    path = fs.join(input_path, arguments.ignore)

    # Inform the user
    log.info("Creating mask covering objects that should be ignored from " + path + " ...")

    # Load the region and create a mask from it
    ignore_region = PixelRegionList.from_file(path)

# No ignore region
else: ignore_region = None

# -----------------------------------------------------------------

# Inform the user
log.info("Running the source finder ...")

# Create a SourceFinder instance
finder = SourceFinder.from_arguments(arguments)

# Run the extractor
finder.run(image.frames.primary, catalog_importer.galactic_catalog, catalog_importer.stellar_catalog, special_region, ignore_region, bad_mask)

# -----------------------------------------------------------------

# Save the galaxy region
galaxy_region = finder.galaxy_region
galaxy_region_path = fs.join(output_path, "galaxies.reg")
galaxy_region.save(galaxy_region_path)

# Save the star region
star_region = finder.star_region
star_region_path = fs.join(output_path, "stars.reg")
star_region.save(star_region_path)

# Save the saturation region
saturation_region = finder.saturation_region
saturation_region_path = fs.join(output_path, "saturation.reg")
saturation_region.save(saturation_region_path)

# Save the other sources region
other_region = finder.other_region
other_region_path = fs.join(output_path, "other_sources.reg")
other_region.save(other_region_path)

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
path = fs.join(output_path, "segments.fits")
segments.save(path)

# -----------------------------------------------------------------

# Get the segmentation maps (galaxies, stars and other sources)
galaxy_segments = finder.galaxy_segments
star_segments = finder.star_segments
other_segments = finder.other_segments

# Get the regions
if arguments.interactive:

    # Wait for a key stroke to continue with the actual extraction
    name = raw_input("Press enter to continue with the extraction step ...")

    # Import the star and saturation regions which have been adjusted by the user
    star_region = PixelRegionList.from_file(star_region_path)
    saturation_region = PixelRegionList.from_file(saturation_region_path)
    other_region = PixelRegionList.from_file(other_region_path)

# -----------------------------------------------------------------

# Inform the user
log.info("Running the source extractor ...")

# Create a SourceExtractor instance
extractor = SourceExtractor.from_arguments(arguments)

# Run the source extractor
extractor.run(image.frames.primary, star_region, saturation_region, galaxy_segments, star_segments, other_segments)

# -----------------------------------------------------------------

# Determine the path to the result
result_path = fs.join(output_path, image.name + ".fits")

# Save the resulting image as a FITS file
image.save(result_path)

# -----------------------------------------------------------------
