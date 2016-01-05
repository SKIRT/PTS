#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.tools.catalogs Contains convenient functions for obtaining information from star or
#  galaxy catalogs.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import astronomical modules
import astropy.units as u
import astropy.coordinates as coord
from astroquery.vizier import Vizier
from astroquery.ned import Ned
from astroquery.irsa_dust import IrsaDust

# Import the relevant AstroMagic classes and modules
from . import regions

# -----------------------------------------------------------------

def galaxies_in_box(center, ra_span, dec_span):

    """
    This function ...
    :return:
    """

    # Initialize a list to contain the galaxies
    names = []

    # Create a new Vizier object and set the row limit to -1 (unlimited)
    viz = Vizier(keywords=["galaxies", "optical"])
    viz.ROW_LIMIT = -1

    # Query Vizier and obtain the resulting table
    result = viz.query_region(center, width=ra_span, height=dec_span, catalog=["VII/237"])

    # I noticed something strange happening once; where there were no entries in the result,
    # with the following parameters:
    #   center = (149.07614359, 69.24847936)
    #   ra_span = 1.600000128 deg
    #   dec_span = 1.3966667784 deg
    #   catalog = ["VII/237"]
    # When ra_span was only slightly changed (e.g. change the last digit to a '7'), output was normal
    # Thus, it seems that the query goes wrong with specific values of the width (and/or height), in which
    # case changing the value very slightly resolves the problem...
    # I am baffled by this and I see no reasonable explanation.
    if len(result) == 0:

        ra_span *= 1.0 + 1e-5
        result = viz.query_region(center, width=ra_span, height=dec_span, catalog=["VII/237"])

    table = result[0]

    # Loop over the rows in the table
    for entry in table:
        name = "PGC " + str(entry["PGC"])
        coordinate = coord.SkyCoord(ra=entry["_RAJ2000"], dec=entry["_DEJ2000"], unit=(u.deg, u.deg), frame='fk5')
        namepluscoordinate = (name, coordinate)
        names.append(namepluscoordinate)

    # Return the list of galaxies
    return names

# -----------------------------------------------------------------

def fetch_objects_in_box(box, catalog, keywords, radius, limit=None, column_filters=None):

    """
    This function ...
    :param box:
    :param catalog:
    :param keywords:
    :param radius:
    :param limit:
    :param column_filters:
    :return:
    """

    # Define the center coordinate for the box
    coordinate = coord.SkyCoord(ra=box[0], dec=box[1], unit=(u.deg, u.deg), frame='fk5') # frame: icrs, fk5... ?

    # Make a Vizier object
    if column_filters is None:
        viz = Vizier(columns=['_RAJ2000', '_DEJ2000','B-V', 'Vmag', 'Plx'], keywords=keywords)
    else:
        viz = Vizier(columns=['_RAJ2000', '_DEJ2000','B-V', 'Vmag', 'Plx'], column_filters=column_filters, keywords=keywords)

    # No limit on the number of entries
    viz.ROW_LIMIT = limit if limit is not None else -1

    # Query the box of our image frame
    result = viz.query_region(coordinate, width=box[3]*u.deg, height=box[2]*u.deg, catalog=catalog)

    region_string = "# Region file format: DS9 version 3.0\n"
    region_string += "global color=green\n"

    # Result may contain multiple tables (for different catalogs)
    for table in result:

        # For every entry in the table
        for entry in table:

            # Get the right ascension and the declination
            ra = entry[0]
            dec = entry[1]

            # Create a string with the coordinates of the star
            regline = "fk5;circle(%s,%s,%.2f\")\n" % (ra, dec, radius)
            #regline = "image;circle(%s,%s,%s)\n" % (ra, dec, radius)

            # Add the parameters of this star to the region string
            region_string += regline

    # Return the region
    return regions.parse(region_string)

# -----------------------------------------------------------------

def fetch_object_by_name(name, radius):

    """
    This function ...
    :param name:
    :param color:
    :param radius:
    :return:
    """

    # Query the NED database for the object
    table = Ned.query_object(name)

    region_string = "# Region file format: DS9 version 3.0\n"
    region_string += "global color=green\n"

    # For every entry in the table
    for entry in table:

        # Get the right ascension and the declination
        ra = entry[2]
        dec = entry[3]

        #print coordinates.degrees_to_hms(ra=ra, dec=dec)

        # Create a string with the coordinates of the star
        regline = "fk5;circle(%s,%s,%.2f\")\n" % (ra, dec, radius)

        # Add the parameters of this star to the region string
        region_string += regline

    # Return the region
    return regions.parse(region_string)

# -----------------------------------------------------------------

def fetch_galactic_extinction(name, filter_name):

    """
    This function ...
    :param name:
    :param filter_name:
    :return:
    """

    table = IrsaDust.get_extinction_table(name)

    for index, item in enumerate(table["Filter_name"]):

        if item == filter_name: break

    return table["A_SandF"][index]

# -----------------------------------------------------------------
