#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.catalog_to_region Create a region file from a catalog file

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import argparse

# Import the relevant PTS classes and modules
from pts.core.tools import tables
from pts.magic.basics.skygeometry import SkyCoordinate
from pts.magic.core.frame import Frame
from pts.magic.tools import statistics

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument("catalog", type=str, help="the name of the region file")
parser.add_argument("image", type=str, help="the name of the image file for which to create the region")
parser.add_argument("fwhm", type=float, help="the FWHM of the stars (in pixels)")
parser.add_argument("sigma_level", type=float, nargs='?', help="the sigma level", default=3.0)
parser.add_argument("--color", type=str, help="the color", default="blue")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Load the catalog
catalog_name = os.path.splitext(os.path.basename(arguments.catalog))[0]
catalog = tables.from_file(arguments.catalog)

# Open the frame
frame = Frame.from_file(arguments.image)

# -----------------------------------------------------------------

# Determine the path to the region file
path = os.path.join(os.getcwd(), catalog_name + ".reg")

# Create a file
f = open(path, 'w')

# Initialize the region string
print("# Region file format: DS9 version 4.1", file=f)

# Create the list of stars
for i in range(len(catalog)):

    # Get the star properties
    #catalog = self.catalog["Catalog"][i]
    star_id = catalog["Id"][i]
    ra = catalog["Right ascension"][i]
    dec = catalog["Declination"][i]
    #ra_error = catalog["Right ascension error"][i] * u.mas
    #dec_error = catalog["Declination error"][i] * u.mas
    #confidence_level = catalog["Confidence level"][i]

    # Create a sky coordinate for the star position
    position = SkyCoordinate(ra=ra, dec=dec, unit="deg", frame="fk5")

    # To pixel position
    center = position.to_pixel(frame.wcs)

    # Determine the radius
    radius = arguments.fwhm * statistics.fwhm_to_sigma * arguments.sigma_level

    # Determine the text (=star index)
    text = star_id.split("/")[1]

    # Show a circle for the star
    suffix = " # "
    color_suffix = "color = " + arguments.color
    text_suffix = "text = {" + text + "}"
    suffix += color_suffix + " " + text_suffix
    print("image;circle({},{},{})".format(center.x+1, center.y+1, radius) + suffix, file=f)

# Close the file
f.close()

# -----------------------------------------------------------------
