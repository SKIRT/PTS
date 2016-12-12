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

# Import astronomical modules
from astropy.units import Unit

# Import the relevant PTS classes and modules
from pts.magic.misc.imageimporter import ImageImporter
from pts.magic.sources.finder import SourceFinder
from pts.magic.catalog.importer import CatalogImporter
from pts.magic.core.image import Image
from pts.magic.region.list import PixelRegionList, SkyRegionList
from pts.core.tools import configuration
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()

# Configuration
parser.add_argument('--config', type=str, help='the name of a configuration file', default=None)
parser.add_argument("--settings", type=configuration.from_string, help="settings")

# Basic
parser.add_argument("image", type=str, help="the name of the input image")

# Options for logging
parser.add_argument("--debug", action="store_true", help="enable debug logging mode")
parser.add_argument('--report', action='store_true', help='write a report file')

# Input and output
parser.add_argument("-i", "--input", type=str, help="the name of the input directory")
parser.add_argument("-o", "--output", type=str, help="the name of the output directory")

# Writing
parser.add_argument("--catalogs", action="store_true", help="save the catalog files")

# Advanced options
parser.add_argument("--principal_region", type=str, help="the path to a region file with a contour of the principal galaxy (in sky coordinates!)")
parser.add_argument("--synchronize", action="store_true", help="synchronize with DustPedia catalog")
parser.add_argument("--filecatalog", action="store_true", help="use file catalogs")
parser.add_argument("--interpolation_method", type=str, help="the interpolation method to use")
parser.add_argument("--downsample", type=float, help="specify the degree of downsampling (no downsampling if not specified)")
parser.add_argument("--no_saturation", action="store_true", help="don't look for saturated stars")
parser.add_argument("--no_other", action="store_true", help="don't look for sources outside of the catalogs")
parser.add_argument("--saturation_dilation_factor", type=float, help="the dilation factor to be used for the detected saturation")
parser.add_argument("--other_dilation_factor", type=float, help="the dilation factor to be used for the detected other sources")
parser.add_argument("--saturation_sigma_level", type=float, help="the sigma level to be used for the segmentation step in saturation detection")
parser.add_argument("--other_sigma_level", type=float, help="the sigma level to be used for the segmentation step in finding the other sources")

parser.add_argument("--stars_peak_sigma_level", type=float, help="the sigma level for the peak finding step of the star finder")

parser.add_argument("--fwhm", type=float, help="the FWHM of the image in arcseconds")

parser.add_argument("--saturation_box_sigmas", type=float, help="ask Sam")

# Input regions
parser.add_argument("--ignore", type=str, help="the name of the file specifying regions to ignore (in sky coordinates!)")
parser.add_argument("--special", type=str, help="the name of the file specifying regions with objects needing special attention (in sky coordinates!)")
parser.add_argument("--bad", type=str, help="the name of the file specifying regions that have to be added to the mask of bad pixels")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# If an input directory is given
if arguments.input is not None:

    # Determine the full path to the input directory
    input_path = fs.absolute_path(arguments.input)

    # Give an error if the input directory does not exist
    if not fs.is_directory(input_path):
        raise argparse.ArgumentError(input_path, "The input directory does not exist")

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
log.start("Starting find_sources ...")

# -----------------------------------------------------------------

# Determine the full path to the image
image_path = fs.absolute_path(arguments.image)

# Determine the full path to the bad region file
bad_region_path = fs.join(input_path, arguments.bad) if arguments.bad is not None else None

# Set FWHM
fwhm = arguments.fwhm * Unit("arcsec") if arguments.fwhm is not None else None

# Import the image
importer = ImageImporter()
importer.run(image_path, bad_region_path=bad_region_path, fwhm=fwhm)

# Get the image
image = importer.image

# Get the mask of bad pixels
bad_mask = image.masks.bad if "bad" in image.masks else None

# -----------------------------------------------------------------

# Create a CatalogImporter instance
catalog_importer = CatalogImporter()

# Set file catalog options
if arguments.filecatalog:

    catalog_importer.config.galaxies.use_catalog_file = True
    catalog_importer.config.galaxies.catalog_path = fs.join(input_path, "galaxies.cat")
    catalog_importer.config.stars.use_catalog_file = True
    catalog_importer.config.stars.catalog_path = fs.join(input_path, "stars.cat")

# Run the catalog importer
catalog_importer.run(image.frames.primary) # work with coordinate box instead ? image.coordinate_box ?

# Save the catalogs if requested
if arguments.catalogs:

    # Determine the full path to the galactic catalog file
    path = fs.join(output_path, "galaxies.cat")

    # Save the galactic catalog
    catalog_importer.write_galactic_catalog_to(path)

    # Determine the full path to the stellar catalog file
    path = fs.join(output_path, "stars.cat")

    # Save the stellar catalog
    catalog_importer.write_stellar_catalog_to(path)

# -----------------------------------------------------------------

# If a special region is defined
if arguments.special is not None:

    # Determine the full path to the special region file
    path = fs.join(input_path, arguments.special)

    # Inform the user
    log.info("Loading region indicating areas that require special attention from " + path + " ...")

    # Load the region and create a mask from it
    special_region = SkyRegionList.from_file(path)

# No special region
else: special_region = None

# -----------------------------------------------------------------

# If an ignore region is defined
if arguments.ignore is not None:

    # Determine the full path to the ignore region file
    path = fs.join(input_path, arguments.ignore)

    # Inform the user
    log.info("Loading region indicating areas that should be ignored from " + path + " ...")

    # Load the region and create a mask from it
    ignore_region = SkyRegionList.from_file(path)

# No ignore region
else: ignore_region = None

# -----------------------------------------------------------------

if arguments.principal_region is not None: arguments.principal_region = fs.join(input_path, arguments.principal_region)

# -----------------------------------------------------------------

# Create a SourceFinder instance
finder = SourceFinder.from_arguments(arguments)

# Run the source finder
finder.run(image.frames.primary, catalog_importer.galactic_catalog, catalog_importer.stellar_catalog, special_region, ignore_region, bad_mask)

# -----------------------------------------------------------------

# Show the FWHM
log.info("The FWHM that could be fitted to the point sources is " + str(finder.fwhm))

# -----------------------------------------------------------------

# Save the galaxy region
galaxy_sky_region = finder.galaxy_sky_region
if galaxy_sky_region is not None:

    galaxy_region = galaxy_sky_region.to_pixel(image.wcs)
    path = fs.join(output_path, "galaxies.reg")
    galaxy_region.save(path)

# Save the star region
star_sky_region = finder.star_sky_region
if star_sky_region is not None:

    star_region = star_sky_region.to_pixel(image.wcs)
    path = fs.join(output_path, "stars.reg")
    star_region.save(path)

# Save the saturation region
saturation_sky_region = finder.saturation_sky_region
if saturation_sky_region is not None:

    saturation_region = saturation_sky_region.to_pixel(image.wcs)
    path = fs.join(output_path, "saturation.reg")
    saturation_region.save(path)

# Save the region of other sources
other_sky_region = finder.other_sky_region
if other_sky_region is not None:

    other_region = other_sky_region.to_pixel(image.wcs)
    path = fs.join(output_path, "other_sources.reg")
    other_region.save(path)

# -----------------------------------------------------------------

# Create an image with the segmentation maps
segments = Image("segments")

# Add the segmentation map of the galaxies
segments.add_frame(finder.galaxy_segments, "galaxies")

# Add the segmentation map of the saturated stars
segments.add_frame(finder.star_segments, "stars")

# Add the segmentation map of the other sources
segments.add_frame(finder.other_segments, "other_sources")

# Save the FITS file with the segmentation maps
path = fs.join(output_path, "segments.fits")
segments.save(path)

# -----------------------------------------------------------------

# Write statistics file
statistics_path = fs.join(output_path, "statistics.dat")
finder.write_statistics(statistics_path)

# -----------------------------------------------------------------
