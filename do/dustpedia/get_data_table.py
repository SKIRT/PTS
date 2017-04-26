#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from pts.dustpedia.core.database import DustPediaDatabase, get_account
from pts.dustpedia.core.sample import DustPediaSample
from pts.core.tools import filesystem as fs
from pts.core.tools.logging import log
from pts.core.basics.table import SmartTable
from pts.core.tools import introspection
from pts.core.tools import tables
from pts.core.filter.filter import parse_filter
from pts.core.units.parsing import parse_unit as u
from pts.core.basics.map import Map

# -----------------------------------------------------------------

path = fs.cwd()

# -----------------------------------------------------------------

calibration_magnitude = dict()
calibration_magnitude["GALEX FUV"] = 0.05
calibration_magnitude["GALEX NUV"] = 0.03
calibration_magnitude["2MASS J"] = 0.03
calibration_magnitude["2MASS H"] = 0.03
calibration_magnitude["2MASS Ks"] = 0.03

# -----------------------------------------------------------------

def calibration_magnitude_for_filter(fltr):

    """
    This function ...
    :param fltr: 
    :return: 
    """

    for key in calibration_magnitude:
        key_filter = parse_filter(key)
        if key_filter == fltr: return calibration_magnitude[key] * u("mag")
    return None

# -----------------------------------------------------------------

# Create the DustPedia sample
sample = DustPediaSample()
galaxy_name = sample.get_name("M81")

# Create the database
database = DustPediaDatabase()

# Login
username, password = get_account()
database.login(username, password)

# -----------------------------------------------------------------

# Create table
table = SmartTable()
table.add_column_info("Observatory/telescope", str, None, "Observatory of telescope")
table.add_column_info("Instrument", str, None, "Instrument")
table.add_column_info("Band", str, None, "Band")
table.add_column_info("Wavelength", float, "micron", "Filter effective wavelength")
table.add_column_info("Pixelscale", float, "arcsec", "Image pixelscale")
table.add_column_info("FWHM of PSF", float, "arcsec", "FWHM of the PSF")
table.add_column_info("Calibration uncertainty (%)", float, None, "Percentual calibration uncertainty")
table.add_column_info("Calibration uncertainty", float, "mag", "Calibration uncertaintity in magnitudes")

# -----------------------------------------------------------------

# Load properties
properties_path = fs.join(introspection.pts_dat_dir("dustpedia"), "properties.dat")
properties = tables.from_file(properties_path, format="ascii.commented_header")

data = OrderedDict()

# Convert properties in dictionary per filter
for index in range(len(properties)):

    facility = properties["Facility"][index]
    band = properties["Band"][index]

    # Skip DSS
    if facility == "DSS1" or facility == "DSS2": continue

    filterstring = facility + " " + band
    fltr = parse_filter(filterstring)

    # Get the properties
    pixelscale = properties["Pixelsize"][index] * u("arcsec")
    fwhm = properties["FWHM"][index] * u("arcsec")
    calibration = properties["Calibration error"][index]

    # Set properties
    data[fltr] = Map(pixelscale=pixelscale, fwhm=fwhm, calibration=calibration)

# -----------------------------------------------------------------

#print([str(key) for key in data.keys()])

# Loop over the images for M81
filters = database.get_image_names_and_filters(galaxy_name)
for name in filters:

    # Get the filter
    fltr = filters[name]

    # Get filter properties
    observatory = fltr.observatory
    instrument = fltr.instrument
    band = fltr.band
    wavelength = fltr.wavelength

    if observatory == "SDSS": observatory = "APO"

    # Get properties
    pixelscale = data[fltr].pixelscale
    fwhm = data[fltr].fwhm
    percentual_calibration = data[fltr].calibration

    # Get magnitude calibration
    calibration = calibration_magnitude_for_filter(fltr)

    # Set row values
    values = [observatory, instrument, band, wavelength, pixelscale, fwhm, percentual_calibration, calibration]

    # Add entry to the table
    table.add_row(values)

# -----------------------------------------------------------------

# Print
table.sort("Wavelength")
#print(table)

# -----------------------------------------------------------------

# Print for latex

header = " & ".join(table.colnames) + " \\\\"
print(header)

units = []
for name in table.colnames:
    unit = table.column_unit(name)
    if unit is None: units.append("")
    else: units.append(str(unit))
units_string = " & ".join(units) + " \\\\"
print(units_string)

for index in range(len(table)):

    row = []
    for name in table.colnames: row.append(str(table[name][index]))
    row_string = " & ".join(row) + " \\\\"

    print(row_string)

# -----------------------------------------------------------------
