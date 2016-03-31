#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.initialize_data Initialize the data for the radiative transfer modeling pipeline.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import argparse

# Import the relevant PTS classes and modules
from pts.magic.core.image import Image
from pts.core.tools import logging, time, filesystem, tables

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()

# Basic options
parser.add_argument("path", type=str, nargs='?', help="the modeling path")

# Logging
parser.add_argument("--debug", action="store_true", help="enable debug logging mode")
parser.add_argument("--report", action='store_true', help='write a report file')

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Set the modeling path
if arguments.path is None: arguments.path = filesystem.cwd()

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = filesystem.join(arguments.path, time.unique_name("log") + ".txt") if arguments.report else None

# Determine the log level
level = "DEBUG" if arguments.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting inspect_data ...")

# -----------------------------------------------------------------

# Determine the path to the data directory
modeling_path = arguments.path
data_path = filesystem.join(modeling_path, "data")

# Initialize the columns of the table
names_column = []
paths_column = []
observatory_column = []
instrument_column = []
band_column = []
unit_column = []
pixelscale_column = []
has_errors_column = []
prep_names_column = []
names = ["Image name", "Image path", "Observatory", "Instrument", "Band", "Unit", "Pixelscale", "Has errors", "Preparation name"]

# Loop over all subdirectories of the data directory
for path, name in filesystem.directories_in_path(data_path, not_contains="bad", returns=["path", "name"]):

    # Get the name of the observatory
    #observatory = name

    # Loop over all FITS files found in the current subdirectory
    for image_path, image_name in filesystem.files_in_path(path, extension="fits", not_contains="_Error", returns=["path", "name"]):

        # Open the image
        image = Image.from_file(image_path)

        # Check if an error map is present (as one of the image frames or as a seperate FITS file)
        has_errors = "errors" in image.frames or filesystem.is_file(filesystem.join(path, image_name + "_Error.fits"))

        if image.filter is not None:
            observatory = image.filter.observatory
            instrument = image.filter.instrument
            band = image.filter.band
            prep_name = instrument + " " + band
        else:
            observatory = None
            instrument = None
            band = None
            prep_name = image_name

        names_column.append(image_name)
        paths_column.append(image_path)
        observatory_column.append(observatory)
        instrument_column.append(instrument)
        band_column.append(band)
        unit_column.append(str(image.unit))
        pixelscale_column.append(image.xy_average_pixelscale.to("arcsec/pix").value)
        has_errors_column.append(has_errors)
        prep_names_column.append(prep_name)

# Create the table
data = [names_column, paths_column, observatory_column, instrument_column, band_column, unit_column, pixelscale_column, has_errors_column, prep_names_column]
table = tables.new(data, names)

# Save the info table
info_path = filesystem.join(data_path, "info.dat")
tables.write(table, info_path)

# -----------------------------------------------------------------
