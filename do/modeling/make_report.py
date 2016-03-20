#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.make_report Make a report of a certain step of the radiative transfer modeling procedure.

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
parser.add_argument("step", type=str, help="the modeling step for which to create the report")
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
logging.setup_log(level=level, path=logfile_path)
logging.log.info("Starting make_report ...")

# -----------------------------------------------------------------

# Data initialization
if arguments.step == "initialization":

    # Determine the path to the preparation directory
    modeling_path = arguments.path
    prep_path = filesystem.join(modeling_path, "prep")

    # Initialize the columns of the table
    names_column = []
    initialized_column = []
    errors_column = []
    galaxy_column = []
    star_column = []
    saturation_column = []
    other_column = []
    segments_column = []
    fwhm_column = []

    names = ["Image name", "Initialized", "Errors", "Galaxy region", "Star region", "Saturation region", "Other sources region", "Segmentation maps", "FWHM"]

    # Loop over all subdirectories of the preparation directory
    for path, name in filesystem.directories_in_path(prep_path, returns=["path", "name"]):

        # Determine the path to the initialized image
        image_path = filesystem.join(path, "initialized.fits")
        initialized = filesystem.is_file(image_path)

        # Determine the path to the sources directory
        sources_path = filesystem.join(path, "sources")
        sources_present = filesystem.is_directory(sources_path)

        # Check the presence of the region files
        has_galaxy_region = filesystem.is_file(filesystem.join(sources_path, "galaxies.reg")) if sources_present else False
        has_star_region = filesystem.is_file(filesystem.join(sources_path, "stars.reg")) if sources_present else False
        has_saturation_region = filesystem.is_file(filesystem.join(sources_path, "saturation.reg")) if sources_present else False
        has_other_region = filesystem.is_file(filesystem.join(sources_path, "other_sources.reg")) if sources_present else False

        # Check the presence of the segmentation image
        has_segments = filesystem.is_file(filesystem.join(sources_path, "segments.fits")) if sources_present else False

        # Open the image
        image = Image.from_file(image_path)

        # Check the presence of an error frame
        has_errors = "errors" in image.frames.keys()

        # Get the FWHM
        fwhm = image.fwhm.to("arcsec").value

        # Fill in the columns
        names_column.append(name)
        initialized_column.append(initialized and sources_present)
        errors_column.append(has_errors)
        galaxy_column.append(has_galaxy_region)
        star_column.append(has_star_region)
        saturation_column.append(has_saturation_region)
        other_column.append(has_other_region)
        segments_column.append(has_segments)
        fwhm_column.append(fwhm)

    # Create the table
    data = [names_column, initialized_column, errors_column, galaxy_column, star_column, saturation_column, other_column, segments_column, fwhm_column]
    table = tables.new(data, names)

    # Save the table
    table_path = filesystem.join(prep_path, "initialization.dat")
    tables.write(table, table_path)

# Data preparation
elif arguments.step == "preparation":

    # Determine the path to the preparation directory
    modeling_path = arguments.path
    prep_path = filesystem.join(modeling_path, "prep")

    # Initialize the columns of the table
    names_column = []
    extracted_column = []
    corrected_column = []
    converted_column = []
    convolved_column = []
    rebinned_column = []
    subtracted_column = []
    result_column = []

    names = ["Image name", "Sources extracted", "Corrected for extinction", "Unit converted", "Convolved", "Rebinned", "Sky subtracted", "Result"]

    # Loop over all subdirectories of the preparation directory
    for path, name in filesystem.directories_in_path(prep_path, returns=["path", "name"]):

        # Determine the path to the extracted image
        extracted_path = filesystem.join(path, "extracted.fits")
        has_extracted = filesystem.is_file(extracted_path)

        # Determine the path to the extinction-corrected image
        corrected_path = filesystem.join(path, "corrected_for_extinction.fits")
        has_corrected = filesystem.is_file(corrected_path)

        # Determine the path to the unit-converted image
        converted_path = filesystem.join(path, "converted_unit.fits")
        has_converted = filesystem.is_file(converted_path)

        # Determine the path to the convolved image
        convolved_path = filesystem.join(path, "convolved.fits")
        has_convolved = filesystem.is_file(convolved_path)

        # Determine the path to the rebinned image
        rebinned_path = filesystem.join(path, "rebinned.fits")
        has_rebinned = filesystem.is_file(rebinned_path)

        # Determine the path to the sky-subtracted image
        subtracted_path = filesystem.join(path, "subtracted.fits")
        has_subtracted = filesystem.is_file(subtracted_path)

        # Determine the path to the prepared image
        result_path = filesystem.join(path, "result.fits")
        has_result = filesystem.is_file(result_path)

        # Fill in the columns
        names_column.append(name)
        extracted_column.append(has_extracted)
        corrected_column.append(has_corrected)
        converted_column.append(has_converted)
        convolved_column.append(has_convolved)
        rebinned_column.append(has_rebinned)
        subtracted_column.append(has_subtracted)
        result_column.append(has_result)

    # Create the table
    data = [names_column, extracted_column, corrected_column, converted_column, convolved_column, rebinned_column, subtracted_column, result_column]
    table = tables.new(data, names)

    # Save the table
    table_path = filesystem.join(prep_path, "preparation.dat")
    tables.write(table, table_path)

# Other
else: raise ValueError("Other modeling steps not enabled yet")

# -----------------------------------------------------------------
