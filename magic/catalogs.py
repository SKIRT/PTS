#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import image modules
import regions
from tools import coordinates

# Import astronomical modules
import astropy.units as u
import astropy.coordinates as coord
from astroquery.vizier import Vizier
from astroquery.ned import Ned

# *****************************************************************

def fetch_objects_in_box(box, catalogs, keywords, radius, limit=None, column_filters=None):

    """
    This function ...
    :param box:
    :param catalogs:
    :param keywords:
    :param radius:
    :param limit:
    :param column_filters:
    :return:
    """

    # Define the center coordinate for the box
    coordinate = coord.SkyCoord(ra=box[0], dec=box[1], unit=(u.deg, u.deg), frame='fk5') # frame: icrs, fk5... ?

    # Make a Vizier object
    viz = Vizier(columns=['_RAJ2000', '_DEJ2000','B-V', 'Vmag', 'Plx'], column_filters=column_filters, keywords=keywords)
    #viz = Vizier(columns=['_RAJ2000', '_DEJ2000'], keywords=keywords)

    # No limit on the number of entries
    viz.ROW_LIMIT = limit if limit is not None else -1

    # Query the box of our image frame
    result = viz.query_region(coordinate, width=box[3]*u.deg, height=box[2]*u.deg, catalog=catalogs)

    region_string = "# Region file format: DS9 version 3.0\n"
    region_string += "global color=green\n"

    # Result may contain multiple tables (for different catalogs)
    for table in result:

        # For every entry in the table
        for entry in table:

            # Get the right ascension and the declination
            ra = entry[0]
            dec = entry[1]

            #print coordinates.degrees_to_hms(ra=ra, dec=dec)

            # Create a string with the coordinates of the star
            regline = "fk5;circle(%s,%s,%.2f\")\n" % (ra, dec, radius)

            # Add the parameters of this star to the region string
            region_string += regline

    # Return the region
    return regions.parse(region_string)

# *****************************************************************

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

# *****************************************************************