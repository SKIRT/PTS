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

# Import standard modules
import copy
import math
import numpy as np

# Import astronomical modules
import astropy.units as u
import astropy.coordinates as coord
from astroquery.vizier import Vizier
from astroquery.ned import Ned
from astroquery.irsa_dust import IrsaDust
import astroquery.exceptions
from astropy.coordinates import Angle

# Import the relevant PTS classes and modules
from ...core.tools import tables

# Import the relevant AstroMagic classes and modules
from . import regions
from ..basics import Position

# -----------------------------------------------------------------

def merge_stellar_catalogs(catalog_a, catalog_b):

    """
    This function ...
    :param catalog_a:
    :param catalog_b:
    :return:
    """

    # Create a copy of catalog a
    new_catalog = copy.deepcopy(catalog_a)

    # Keep track of which stars in catalog a have been encountered in catalog b (to avoid unnecessary checking)
    encountered = [False] * len(catalog_a)

    # Loop over the entries in catalog b
    for j in range(len(catalog_b)):

        # Loop over the entries in catalog a to find a match
        for i in range(len(catalog_a)):

            # Skip this entry (star) if it has already been matched with a star from the other catalog
            if encountered[i]: continue

            # If star i in catalog a is the same as star j in catalog b, break the loop over catalog a's stars
            if catalog_a["Catalog"][i] == catalog_b["Catalog"][j] and catalog_a["Id"][i] == catalog_b["Id"][j]:

                encountered[i] = True
                break

        # If a break is not encountered, a match is not found -> add the star from catalog b
        else:

            # Add the corresponding row
            new_catalog.add_row(catalog_b[j])

    # Return the new catalog
    return new_catalog

# -----------------------------------------------------------------

def merge_galactic_catalogs(catalog_a, catalog_b):

    """
    This function ...
    :param catalog_a:
    :param catalog_b:
    :return:
    """

    # Create a copy of catalog a
    new_catalog = copy.deepcopy(catalog_a)

    # Keep track of which galaxies in catalog a have been encountered in catalog b (to avoid unnecessary checking)
    encountered = [False] * len(catalog_a)

    # Loop over the entries in catalog b
    for j in range(len(catalog_b)):

        # Loop over the entries in catalog a to find a match
        for i in range(len(catalog_a)):

            # Skip this entry (galaxy) if it has already been matches with a galaxy from the other catalog
            if encountered[i]: continue

            # If galaxy i in catalog a is the same as galaxy j in catalog b, break the loop over catalog a's galaxies
            if catalog_a["Name"][i] == catalog_b["Name"][j]:

                encountered[i] = True
                break

        # If a break is not encountered, a match is not found -> add the galaxy from catalog b
        else:

            # Add the corresponding row
            new_catalog.add_row(catalog_b[j])

    # Return the new catalog
    return new_catalog

# -----------------------------------------------------------------

def from_stars(stars):

    """
    This function ...
    :param stars:
    :return:
    """

    # Initialize empty lists for the table columns
    catalog_column = []
    id_column = []
    ra_column = []
    dec_column = []
    ra_error_column = []
    dec_error_column = []
    on_galaxy_column = []
    confidence_level_column = []

    # Loop over all stars
    for star in stars:

        # Fill in the columns with the star properties
        catalog_column.append(star.catalog)
        id_column.append(star.id)
        ra_column.append(star.position.ra.value)
        dec_column.append(star.position.dec.value)
        ra_error_column.append(star.ra_error.value)
        dec_error_column.append(star.dec_error.value)
        on_galaxy_column.append(star.on_galaxy)
        confidence_level_column.append(star.confidence_level)

    # Create and return the table
    data = [catalog_column, id_column, ra_column, dec_column, ra_error_column, dec_error_column, confidence_level_column]
    names = ['Catalog', 'Id', 'Right ascension', 'Declination', 'Right ascension error', 'Declination error', 'Confidence level']

    # Create the catalog
    meta = {'name': 'stars'}
    catalog = tables.new(data, names, meta)

    # Set units
    catalog["Right ascension"].unit = "deg"
    catalog["Declination"].unit = "deg"
    catalog["Right ascension error"].unit = "mas"
    catalog["Declination error"].unit = "mas"

    # Return the catalog
    return catalog

# -----------------------------------------------------------------

def from_galaxies(galaxies):

    """
    This function ...
    :param galaxies:
    :return:
    """

    # Initialize empty lists for the table columns
    name_column = []
    ra_column = []
    dec_column = []
    redshift_column = []
    type_column = []
    alternative_names_column = []
    distance_column = []
    inclination_column = []
    d25_column = []
    major_column = []
    minor_column = []
    pa_column = []

    principal_column = []
    companions_column = []
    parent_column = []

    # Loop over all galaxies
    for galaxy in galaxies:

        # Fill in the columns with the galaxy properties
        pass

    # Create the data structure and names list
    data = [name_column, ra_column, dec_column, redshift_column, type_column, alternative_names_column, distance_column,
            inclination_column, d25_column, major_column, minor_column, pa_column, principal_column, companions_column,
            parent_column]
    names = ["Name", "Right ascension", "Declination", "Redshift", "Type", "Alternative names", "Distance",
             "Inclination", "D25", "Major axis length", "Minor axis length", "Position angle", "Principal",
             "Companion galaxies", "Parent galaxy"]
    meta = {'name': 'stars'}

    # Create the catalog table
    catalog = tables.new(data, names, meta)

    # Set the column units
    catalog["Distance"].unit = "Mpc"
    catalog["Inclination"].unit = "deg"
    catalog["D25"].unit = "arcmin"
    catalog["Major axis length"].unit = "arcmin"
    catalog["Minor axis length"].unit = "arcmin"
    catalog["Position angle"].unit = "deg"
    catalog["Right ascension"].unit = "deg"
    catalog["Declination"].unit = "deg"

    # Return the catalog
    return catalog

# -----------------------------------------------------------------

def create_star_catalog(frame, catalogs=None):

    """
    This function ...
    :return:
    """

    # Initialize empty lists for the table columns
    catalog_column = []
    id_column = []
    ra_column = []
    dec_column = []
    ra_error_column = []
    dec_error_column = []
    magnitude_columns = {}
    magnitude_error_columns = {}
    on_galaxy_column = []
    confidence_level_column = []

    # Get the range of right ascension and declination of this image
    try: center, ra_span, dec_span = frame.coordinate_range()
    except AssertionError as error:

        print("WARNING: The coordinate system and pixelscale do not match")
        center, ra_span, dec_span = frame.coordinate_range(silent=True)

    # Create a new Vizier object and set the row limit to -1 (unlimited)
    viz = Vizier(keywords=["stars", "optical"])
    viz.ROW_LIMIT = -1

    # Loop over the different catalogs
    for catalog in catalogs:

        # Initialize a list to specify which of the stars added to the columns from other catalogs is already
        # matched to a star of the current catalog
        encountered = [False] * len(catalog_column)

        # Inform the user
        print("DEBUG: Querying the " + catalog + " catalog")

        # Query Vizier and obtain the resulting table
        result = viz.query_region(center, width=ra_span, height=dec_span, catalog=catalog)
        table = result[0]

        number_of_stars = 0
        number_of_stars_in_frame = 0
        number_of_new_stars = 0

        magnitudes = {}
        magnitude_errors = {}

        # Get the magnitude in different bands
        for name in table.colnames:

            # If this column name does not end with "mag", skip it
            if not name.endswith("mag"): continue

            # If the column name contains more than one character before "mag", skip it
            if len(name.split("mag")[0]) > 1: continue

            # Get the name of the band
            band = name.split("mag")[0]

            # Create empty lists for the magnitudes and errors
            magnitudes[band] = []
            magnitude_errors[band] = []

        # Loop over all entries in the table
        for i in range(len(table)):

            # -- General information --

            # Get the ID of this star in the catalog
            if catalog == "UCAC4": star_id = table["UCAC4"][i]
            elif catalog == "NOMAD": star_id = table["NOMAD1"][i]
            elif catalog == "II/246": star_id = table["_2MASS"][i]
            else: raise ValueError("Catalogs other than 'UCAC4', 'NOMAD' or 'II/246' are currently not supported")

            # -- Positional information --

            # Get the position of the star as a SkyCoord object and as pixel coordinate
            position = coord.SkyCoord(ra=table["_RAJ2000"][i], dec=table["_DEJ2000"][i], unit=(u.deg, u.deg), frame='fk5')
            pixel_position_x, pixel_position_y = position.to_pixel(frame.wcs, origin=0, mode='wcs')
            pixel_position = Position(pixel_position_x, pixel_position_y)

            # Get the right ascension and declination for the current star
            star_ra = table["_RAJ2000"][i]
            star_dec = table["_DEJ2000"][i]

            number_of_stars += 1

            # If this star does not lie within the frame, skip it
            if not frame.contains(position): continue

            number_of_stars_in_frame += 1

            # Get the mean error on the right ascension and declination
            if catalog == "UCAC4" or catalog == "NOMAD":

                ra_error = table["e_RAJ2000"][i] * u.mas
                dec_error = table["e_DEJ2000"][i] * u.mas

            elif catalog == "II/246":

                error_maj = table["errMaj"][i] * u.arcsec
                error_min = table["errMin"][i] * u.arcsec
                error_theta = Angle(table["errPA"][i], u.deg)

                # Temporary: use only the major axis error (convert the error ellipse into a circle)
                ra_error = error_maj.to("mas")
                dec_error = error_maj.to("mas")

            else: raise ValueError("Catalogs other than 'UCAC4', 'NOMAD' or 'II/246' are currently not supported")

            # -- Magnitudes --

            # Loop over the different bands for which a magnitude is defined
            for band in magnitudes:

                # Determine the column name
                column_name = band + "mag"

                value = table[column_name][i]

                if isinstance(value, np.ma.core.MaskedConstant):

                    magnitudes[band].append(None)
                    magnitude_errors[band].append(None)

                else:

                    # Add the magnitude value
                    magnitudes[band].append(u.Magnitude(value))

                    # Check for presence of error on magnitude
                    error_column_name = "e_" + column_name
                    if error_column_name in table.colnames:
                        error = table[error_column_name][i]
                        if isinstance(error, np.ma.core.MaskedConstant): magnitude_errors[band].append(None)
                        else: magnitude_errors[band].append(u.Magnitude(error))
                    else: magnitude_errors[band].append(None)

            # -- Cross-referencing with previous catalogs --

            # If there are already stars in the list, check for correspondences with the current stars
            for index in range(len(encountered)):

                # Skip stars that are already encountered as matches with the current catalog (we assume there can only
                # be one match of a star of one catalog with the star of another catalog, within the radius of 3 pixels)
                if encountered[index]: continue

                saved_star_position = coord.SkyCoord(ra=ra_column[index], dec=dec_column[index], unit=(u.deg, u.deg), frame='fk5')
                saved_star_pixel_position_x, saved_star_pixel_position_y = saved_star_position.to_pixel(frame.wcs, origin=0, mode="wcs")
                saved_star_pixel_position = Position(saved_star_pixel_position_x, saved_star_pixel_position_y)

                # Calculate the distance between the star already in the list and the new star
                difference = saved_star_pixel_position - pixel_position

                # Check whether the distance is less then 3 pixels
                if difference.norm < 3.0:

                    # Inform the user
                    #print("DEBUG: Star " + star_id + " could be identified with star " + id_column[index] + " from the " + catalog_column[index] + " catalog")

                    # Increment the confidence level for the 'saved' star
                    confidence_level_column[index] += 1

                    # Set the 'encountered' flag to True for the 'saved' star
                    encountered[index] = True

                    # Break, because the current star does not have to be saved again (it is already in the lists)
                    break

            # If no other stars are in the list yet or no corresponding star was found (no break was
            # encountered), just add all stars of the current catalog
            else:

                number_of_new_stars += 1

                # Inform the user
                #print("DEBUG: Adding star " + star_id + " at " + str(position.to_string("hmsdms")))

                # Fill in the column lists
                catalog_column.append(catalog)
                id_column.append(star_id)
                ra_column.append(star_ra)
                dec_column.append(star_dec)
                ra_error_column.append(ra_error.value)
                dec_error_column.append(dec_error.value)
                confidence_level_column.append(1)

        # Debug messages
        print("DEBUG: Number of stars that were in the catalog: " + str(number_of_stars))
        print("DEBUG: Number of stars that fell within the frame: " + str(number_of_stars_in_frame))
        print("DEBUG: Number of stars that were only present in this catalog: " + str(number_of_new_stars))

    # Create and return the table
    data = [catalog_column, id_column, ra_column, dec_column, ra_error_column, dec_error_column, confidence_level_column]
    names = ['Catalog', 'Id', 'Right ascension', 'Declination', 'Right ascension error', 'Declination error', 'Confidence level']

    # TODO: add magnitudes to the table ?

    #magnitude_column_names = []
    #for band in magnitudes:

        # Values
        ##column = MaskedColumn(magnitudes[band], mask=[mag is None for mag in magnitudes[band]])
        ##data.append(column)
        #data.append(magnitudes[band])
        #column_name = band + " magnitude"
        #names.append(column_name)
        #magnitude_column_names.append(column_name)

        # Errors
        ##column = MaskedColumn(magnitude_errors[band], mask=[mag is None for mag in magnitude_errors[band]])
        ##data.append(column)
        #data.append(magnitude_errors[band])
        #column_name = band + " magnitude error"
        #names.append(column_name)
        #magnitude_column_names.append(column_name)

    # Create the catalog
    meta = {'name': 'stars'}
    catalog = tables.new(data, names, meta)

    # Set units
    catalog["Right ascension"].unit = "deg"
    catalog["Declination"].unit = "deg"
    catalog["Right ascension error"].unit = "mas"
    catalog["Declination error"].unit = "mas"
    #for name in magnitude_column_names:
    #    self.catalog[name].unit = "mag"

    # Return the catalog
    return catalog

# -----------------------------------------------------------------

def create_galaxy_catalog(frame):

    """
    This function ...
    :return:
    """

    # Initialize empty lists for the table columns
    name_column = []
    ra_column = []
    dec_column = []
    redshift_column = []
    type_column = []
    alternative_names_column = []
    distance_column = []
    inclination_column = []
    d25_column = []
    major_column = []
    minor_column = []
    pa_column = []

    # Get the range of right ascension and declination of the image
    try: center, ra_span, dec_span = frame.coordinate_range()
    except AssertionError as error:

        print("WARNING: The coordinate system and pixelscale do not match")
        #print(error)
        center, ra_span, dec_span = frame.coordinate_range(silent=True)

    # Find galaxies in the box defined by the center and RA/DEC ranges
    for name, position in galaxies_in_box(center, ra_span, dec_span):

        # Get galaxy information
        gal_name, position, gal_redshift, gal_type, gal_names, gal_distance, gal_inclination, gal_d25, gal_major, gal_minor, gal_pa = get_galaxy_info(name, position)

        # Calculate pixel position in the frame
        pixel_position_x, pixel_position_y = position.to_pixel(frame.wcs, origin=0, mode='wcs')
        pixel_position = Position(pixel_position_x, pixel_position_y)

        # Check whether the pixel position falls within the frame
        if pixel_position.x < 0.0 or pixel_position.x >= frame.xsize: continue
        if pixel_position.y < 0.0 or pixel_position.y >= frame.ysize: continue

        # Fill the columns
        name_column.append(gal_name)
        ra_column.append(position.ra.value)
        dec_column.append(position.dec.value)
        redshift_column.append(gal_redshift)
        type_column.append(gal_type)
        alternative_names_column.append(", ".join(gal_names) if len(gal_names) > 0 else None)
        distance_column.append(gal_distance.value if gal_distance is not None else None)
        inclination_column.append(gal_inclination.degree if gal_inclination is not None else None)
        d25_column.append(gal_d25.value if gal_d25 is not None else None)
        major_column.append(gal_major.value if gal_major is not None else None)
        minor_column.append(gal_minor.value if gal_minor is not None else None)
        pa_column.append(gal_pa.degree if gal_pa is not None else None)

    # Determine the number of galaxies in the lists
    number_of_galaxies = len(name_column)

    # Indicate which galaxy is the principal galaxy
    principal_column = [False] * number_of_galaxies
    principal_index = max(range(number_of_galaxies), key=lambda index: major_column[index])
    principal_column[principal_index] = True

    # Loop over the other galaxies, check if they are companion galaxies of the principal galax
    companions_list_list = [[] for _ in range(number_of_galaxies)]
    parent_column = [None for _ in range(number_of_galaxies)]
    for j in range(number_of_galaxies):

        # Skip the principal galaxy
        if j == principal_index: continue

        # By default, set the 'companion' flag to False
        companion = False

        # Set HII regions as companion galaxies
        if type_column[j] == "HII": companion = True

        # If the galaxies differ in name because of an 'a', 'b' at the end
        if name_column[j][:-1].lower() == name_column[principal_index][:-1].lower(): companion = True

        # If the current galaxy is a companion galaxy
        if companion:

            companions_list_list[principal_index].append(name_column[j])
            parent_column[j] = name_column[principal_index]

    companions_column = []
    for companion_list in companions_list_list:
        if len(companion_list) == 0: companions_column.append(None)
        else: companions_column.append(", ".join(companion_list))

    # Create the data structure and names list
    data = [name_column, ra_column, dec_column, redshift_column, type_column, alternative_names_column, distance_column,
            inclination_column, d25_column, major_column, minor_column, pa_column, principal_column, companions_column,
            parent_column]
    names = ["Name", "Right ascension", "Declination", "Redshift", "Type", "Alternative names", "Distance",
             "Inclination", "D25", "Major axis length", "Minor axis length", "Position angle", "Principal",
             "Companion galaxies", "Parent galaxy"]
    meta = {'name': 'stars'}

    # Create the catalog table
    catalog = tables.new(data, names, meta)

    # Set the column units
    catalog["Distance"].unit = "Mpc"
    catalog["Inclination"].unit = "deg"
    catalog["D25"].unit = "arcmin"
    catalog["Major axis length"].unit = "arcmin"
    catalog["Minor axis length"].unit = "arcmin"
    catalog["Position angle"].unit = "deg"
    catalog["Right ascension"].unit = "deg"
    catalog["Declination"].unit = "deg"

    # Return the catalog
    return catalog

# -----------------------------------------------------------------

def get_galaxy_info(name, position):

    """
    This function ...
    :param name:
    :param position:
    :return:
    """

    # Obtain more information about this galaxy
    try:

        ned_result = Ned.query_object(name)
        ned_entry = ned_result[0]

        # Get a more common name for this galaxy (sometimes, the name obtained from NED is one starting with 2MASX .., use the PGC name in this case)
        if ned_entry["Object Name"].startswith("2MASX "): gal_name = name
        else: gal_name = ned_entry["Object Name"]

        # Get the redshift
        gal_redshift = ned_entry["Redshift"]
        if isinstance(gal_redshift, np.ma.core.MaskedConstant): gal_redshift = None

        # Get the type (G=galaxy, HII ...)
        gal_type = ned_entry["Type"]
        if isinstance(gal_type, np.ma.core.MaskedConstant): gal_type = None

    except astroquery.exceptions.RemoteServiceError:

        # Set attributes
        gal_name = name
        gal_redshift = None
        gal_type = None

    # Create a new Vizier object and set the row limit to -1 (unlimited)
    viz = Vizier(keywords=["galaxies", "optical"])
    viz.ROW_LIMIT = -1

    # Query Vizier and obtain the resulting table
    result = viz.query_object(name.replace(" ", ""), catalog=["VII/237"])

    # Not found ... TODO: fix this ... this object was in the first query output
    if len(result) == 0: return name, position, None, None, [], None, None, None, None, None, None

    table = result[0]

    # Get the correct entry (sometimes, for example for mergers, querying with the name of one galaxy gives two hits! We have to obtain the right one each time!)
    if len(table) == 0: raise ValueError("The galaxy could not be found under this name")
    elif len(table) == 1: entry = table[0]
    else:

        entry = None

        # Some rows don't have names, if no match is found based on the name just take the row that has other names defined
        rows_with_names = []
        for row in table:
            if row["ANames"]: rows_with_names.append(row)

        # If only one row remains, take that one for the galaxy we are looking for
        if len(rows_with_names) == 1: entry = row

        # Else, loop over the rows where names are defined and look for a match
        else:
            for row in rows_with_names:

                names = row["ANames"]

                if name.replace(" ", "") in names or gal_name.replace(" ", "") in names:

                    entry = row
                    break

        # If no matches are found, look for the table entry for which the coordinate matches the given position (if any)
        if position is not None:
            for row in table:
                if row["_RAJ2000"] == position.ra.value and row["_DEJ2000"] == position.dec.value:
                    entry = row
                    break

    # Get the right ascension and the declination
    position = coord.SkyCoord(ra=entry["_RAJ2000"], dec=entry["_DEJ2000"], unit=(u.deg, u.deg), frame='fk5')

    # Get the names given to this galaxy
    gal_names = entry["ANames"].split() if entry["ANames"] else []

    # Get the size of the galaxy
    ratio = np.power(10.0, entry["logR25"]) if entry["logR25"] else None
    diameter = np.power(10.0, entry["logD25"]) * 0.1 * u.arcmin if entry["logD25"] else None

    #print("  D25_diameter = ", diameter)

    radial_profiles_result = viz.query_object(name, catalog="J/ApJ/658/1006")

    if len(radial_profiles_result) > 0:

        radial_profiles_entry = radial_profiles_result[0][0]

        gal_distance = radial_profiles_entry["Dist"] * u.Unit("Mpc")
        gal_inclination = Angle(radial_profiles_entry["i"], u.deg)
        gal_d25 = radial_profiles_entry["D25"] * u.arcmin

    else:

        gal_distance = None
        gal_inclination = None
        gal_d25 = None

    # Get the size of major and minor axes
    gal_major = diameter
    gal_minor = diameter / ratio if diameter is not None and ratio is not None else None

    # Get the position angle of the galaxy
    gal_pa = Angle(entry["PA"] - 90.0, u.deg) if entry["PA"] else None

    # Create and return a new Galaxy instance
    return gal_name, position, gal_redshift, gal_type, gal_names, gal_distance, gal_inclination, gal_d25, gal_major, gal_minor, gal_pa

# -----------------------------------------------------------------

def galaxies_in_box(center, ra_span, dec_span):

    """
    This function ...
    :return:
    """

    # Initialize a list to contain the galaxies
    names = []

    # Other way ?? Much more results ?
    #ra_radius = 0.5 * ra_span.value
    #dec_radius = 0.5 * dec_span.value
    #radius = math.sqrt(ra_radius**2 + dec_radius**2)
    #result_table = Ned.query_region(center, radius=radius)

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
