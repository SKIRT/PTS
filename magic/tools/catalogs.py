#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.tools.catalogs Contains convenient functions for obtaining information from star or
#  galaxy catalogs.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import copy
import requests
import numpy as np
from lxml import html

# Import astronomical modules
from astropy.units import Magnitude
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
from astroquery.ned import Ned
import astroquery.exceptions
from astropy.coordinates import Angle
from astropy.units import UnrecognizedUnit

# Import the relevant PTS classes and modules
from ...core.tools import tables
from ...core.basics.log import log
from ..basics.coordinate import SkyCoordinate
from ..basics.vector import Extent
from ...core.units.parsing import parse_unit as u
from ...core.tools.stringify import tostr

# -----------------------------------------------------------------

leda_search_object_url = "http://leda.univ-lyon1.fr/ledacat.cgi?"

# -----------------------------------------------------------------

stellar_catalogs = dict()
stellar_catalogs["GSC2.3"] = "I/305/out"
stellar_catalogs["NOMAD"] = "I/297/out"
stellar_catalogs["PPMX"] = "I/312"
stellar_catalogs["SDSS9"] = "V/139"
stellar_catalogs["UCAC4"] = "I/322A/out"
stellar_catalogs["URAT1"] = "I/329"
stellar_catalogs["2MASS"] = "II/246/out"

# -----------------------------------------------------------------

stellar_catalog_descriptions = dict()
stellar_catalog_descriptions["GSC2.3"] = "The Guide Star Catalog, Version 2.3.2"
stellar_catalog_descriptions["NOMAD"] = "NOMAD-1 Catalog"
stellar_catalog_descriptions["PPMX"] = "PPMX Catalog of positions and proper motions"
stellar_catalog_descriptions["SDSS9"] = "SDSS Release 9"
stellar_catalog_descriptions["UCAC4"] = "Fourth U.S. Naval Observatory CCD Astrograph Catalog"
stellar_catalog_descriptions["URAT1"] = "URAT1 catalog"
stellar_catalog_descriptions["2MASS"] = "2MASS Point Sources"

# -----------------------------------------------------------------

stellar_catalog_star_id_colname = dict()
stellar_catalog_star_id_colname["GSC2.3"] = "GSC2.3"
stellar_catalog_star_id_colname["NOMAD"] = "NOMAD1"
stellar_catalog_star_id_colname["PPMX"] = "PPMX"
stellar_catalog_star_id_colname["SDSS9"] = "SDSS9"
stellar_catalog_star_id_colname["UCAC4"] = "UCAC4"
stellar_catalog_star_id_colname["URAT1"] = "URAT1"
stellar_catalog_star_id_colname["2MASS"] = "_2MASS"

# -----------------------------------------------------------------

def get_star_id(catalog_name, catalog, index):

    """
    This function ...
    :param catalog_name:
    :param catalog:
    :param index:
    :return:
    """

    colname = get_star_id_colname_for_name(catalog_name)
    return catalog[colname][index]

# -----------------------------------------------------------------

def get_star_id_colname_for_name(name):

    """
    This function ...
    :param name:
    :return:
    """

    return stellar_catalog_star_id_colname[name]

# -----------------------------------------------------------------

def get_stellar_catalog_names():

    """
    This function ...
    :return:
    """

    return stellar_catalogs.keys()

# -----------------------------------------------------------------

def get_stellar_catalog_codes():

    """
    This function ...
    :return:
    """

    return stellar_catalogs.values()

# -----------------------------------------------------------------

def get_stellar_catalog_code_for_name(name):

    """
    This function ...
    :param name:
    :return:
    """

    return stellar_catalogs[name]

# -----------------------------------------------------------------

def get_stellar_catalog_name_for_code(code):

    """
    This function ...
    :param code:
    :return:
    """

    return stellar_catalogs.keys()[stellar_catalogs.values().index(code)]

# -----------------------------------------------------------------

def get_hyperleda_name(galaxy_name):

    """
    This function ...
    :param galaxy_name:
    :return:
    """

    url = leda_search_object_url + galaxy_name

    page_as_string = requests.get(url).content

    tree = html.fromstring(page_as_string)

    tables = [e for e in tree.iter() if e.tag == 'table']

    table = tables[1]
    #table = tables[1]

    #print(table.text_content())
    objname = table.text_content().split(" (")[0]

    #table_rows = [e for e in table.iter() if e.tag == 'tr']
    #column_headings = [e.text_content() for e in table_rows[0].iter() if e.tag == 'th']

    # return table_rows, column_headings

    #objname = str(table_rows[0].text_content().split("\n")[1]).strip()

    # Return the HYPERLEDA name
    return objname

# -----------------------------------------------------------------

def get_ngc_name(galaxy_name, delimiter=" "):

    """
    This function ...
    :param galaxy_name:
    :param delimiter:
    :return:
    """

    # The Simbad querying object
    simbad = Simbad()
    simbad.ROW_LIMIT = -1

    result = simbad.query_objectids(galaxy_name)
    if result is None or len(result) == 0: raise ValueError("Galaxy name '" + galaxy_name + "' could not be recognized")

    # Loop over the results
    for name in result["ID"]:

        if "NGC" in name:

            splitted = name.split("NGC")
            if splitted[0] == "":
                number = int(splitted[1])
                return "NGC" + delimiter + str(number)

    # If nothing is found, return None
    return None

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

            exact_catalog_match = catalog_a["Catalog"][i] == catalog_b["Catalog"][j] and catalog_a["Id"][i] == catalog_b["Id"][j]
            original_catalog_match = catalog_a["Original catalog and id"] == catalog_b["Original catalog and id"]

            # If star i in catalog a is the same as star j in catalog b, break the loop over catalog a's stars
            if exact_catalog_match or original_catalog_match:

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
    meta = {'name': 'galaxies'}

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

def create_star_catalog(coordinate_box, pixelscale, catalogs, check_in_box=False):

    """
    This function ...
    :param coordinate_box:
    :param pixelscale:
    :param catalogs:
    :param check_in_box:
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
    #center, ra_span, dec_span = frame.coordinate_range

    center = coordinate_box.center
    ra_span = 2.0 * coordinate_box.radius.ra
    dec_span = 2.0 * coordinate_box.radius.dec

    # Create a new Vizier object and set the row limit to -1 (unlimited)
    viz = Vizier(keywords=["stars", "optical"])
    viz.ROW_LIMIT = -1

    # Loop over the different catalogs
    for catalog in catalogs:

        # Get catalog code
        code = get_stellar_catalog_code_for_name(catalog)

        # Initialize a list to specify which of the stars added to the columns from other catalogs is already
        # matched to a star of the current catalog
        encountered = [False] * len(catalog_column)

        # Inform the user
        log.debug("Querying the " + catalog + " catalog ...")

        # Query Vizier and obtain the resulting table
        result = viz.query_region(center.to_astropy(), width=ra_span, height=dec_span, catalog=code)
        if len(result) == 0:
            log.warning("No point sources could be found around " + str(center) + " with a RA span of " + str(ra_span) + " and a DEC span of " + str(dec_span) + " with the '" + catalog + "' catalog")
            continue
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

            # Debugging
            log.debug("Processing entry " + str(i+1) + " ...")

            # -- General information --

            # Get the ID of the star in the catalog
            star_id = get_star_id(catalog, table, i)

            # -- Positional information --

            # Get the position of the star as a SkyCoord object and as pixel coordinate
            position = SkyCoordinate(ra=table["RAJ2000"][i], dec=table["DEJ2000"][i], unit="deg", frame="fk5")
            #pixel_position = position.to_pixel(frame.wcs)

            # Get the right ascension and declination for the current star
            star_ra = table["RAJ2000"][i]
            star_dec = table["DEJ2000"][i]

            number_of_stars += 1

            # Optional because this takes a lot of time
            if check_in_box:

                # If this star does not lie within the frame, skip it
                #if not frame.contains(position): continue
                if not coordinate_box.contains(position): continue

            number_of_stars_in_frame += 1

            # DOESN'T WORK ANYMORE!
            # Get the mean error on the right ascension and declination
            #if catalog == "UCAC4" or catalog == "NOMAD":
            #    ra_error = table["e_RAJ2000"][i] * u("mas")
            #    dec_error = table["e_DEJ2000"][i] * u("mas")
            #elif catalog == "II/246":
            #    error_maj = table["errMaj"][i] * u("arcsec")
            ##    error_min = table["errMin"][i] * u("arcsec")
            #    error_theta = Angle(table["errPA"][i], "deg")
            #    # Temporary: use only the major axis error (convert the error ellipse into a circle)
            #    ra_error = error_maj.to("mas")
            #    dec_error = error_maj.to("mas")
            #else: raise ValueError("Catalogs other than 'UCAC4', 'NOMAD' or 'II/246' are currently not supported")
            ra_error = dec_error = None

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
                    magnitudes[band].append(Magnitude(value))

                    # Check for presence of error on magnitude
                    error_column_name = "e_" + column_name
                    if error_column_name in table.colnames:
                        error = table[error_column_name][i]
                        if isinstance(error, np.ma.core.MaskedConstant): magnitude_errors[band].append(None)
                        else: magnitude_errors[band].append(Magnitude(error))
                    else: magnitude_errors[band].append(None)

            # -- Cross-referencing with previous catalogs --

            # If there are already stars in the list, check for correspondences with the current stars
            for index in range(len(encountered)):

                # Skip stars that are already encountered as matches with the current catalog (we assume there can only
                # be one match of a star of one catalog with the star of another catalog, within the radius of 3 pixels)
                if encountered[index]: continue

                saved_star_position = SkyCoordinate(ra=ra_column[index].value, dec=dec_column[index].value, unit="deg", frame="fk5")
                #saved_star_pixel_position = saved_star_position.to_pixel(frame.wcs)

                # Calculate the distance between the star already in the list and the new star
                #difference = saved_star_pixel_position - pixel_position

                difference_ra = saved_star_position.ra - position.ra
                difference_dec = saved_star_position.dec - position.dec
                difference = Extent((difference_ra / pixelscale.to("deg")).to("").value, (difference_dec / pixelscale.to("deg")).to("").value)

                # Check whether the distance is less then 3 pixels
                if difference.norm < 3.0:

                    # Inform the user
                    log.debug("Star " + star_id + " could be identified with star " + id_column[index] + " from the " + catalog_column[index] + " catalog")

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
                ra_column.append(star_ra * u("deg"))
                dec_column.append(star_dec * u("deg"))
                ra_error_column.append(ra_error.value if ra_error is not None else None)
                dec_error_column.append(dec_error.value if dec_error is not None else None)
                confidence_level_column.append(1)

        # Debug messages
        log.debug("Number of stars that were in the catalog: " + str(number_of_stars))
        log.debug("Number of stars that fell within the frame: " + str(number_of_stars_in_frame))
        log.debug("Number of stars that were only present in this catalog: " + str(number_of_new_stars))

    # Create and return the table
    #data = [catalog_column, id_column, ra_column, dec_column, ra_error_column, dec_error_column, confidence_level_column]
    #names = ['Catalog', 'Id', 'Right ascension', 'Declination', 'Right ascension error', 'Declination error', 'Confidence level']

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
    #meta = {'name': 'stars'}
    #catalog = tables.new(data, names, meta)

    # Set units
    #catalog["Right ascension"].unit = "deg"
    #catalog["Declination"].unit = "deg"
    #catalog["Right ascension error"].unit = "mas"
    #catalog["Declination error"].unit = "mas"
    #for name in magnitude_column_names:
    #    self.catalog[name].unit = "mag"

    # Return the catalog
    #return catalog

    return catalog_column, id_column, ra_column, dec_column, ra_error_column, dec_error_column, confidence_level_column

# -----------------------------------------------------------------

def create_galaxy_catalog(coordinate_box):

    """
    This function ...
    :param coordinate_box:
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
    #center, ra_span, dec_span = frame.coordinate_range

    center = coordinate_box.center
    ra_span = 2.0 * coordinate_box.radius.ra
    dec_span = 2.0 * coordinate_box.radius.dec

    # Find galaxies in the box defined by the center and RA/DEC ranges
    index = 0
    for name, position in galaxies_in_box(center, ra_span, dec_span):

        # Debugging
        log.debug("Processing entry " + str(index+1) + " ...")

        index += 1

        # Get galaxy information
        gal_name, position, gal_redshift, gal_type, gal_names, gal_distance, gal_inclination, gal_d25, gal_major, gal_minor, gal_pa = get_galaxy_info(name, position)

        # Calculate pixel position in the frame
        #pixel_position = position.to_pixel(frame.wcs)

        # Check whether the pixel position falls within the frame
        #if pixel_position.x < 0.0 or pixel_position.x >= frame.xsize: continue
        #if pixel_position.y < 0.0 or pixel_position.y >= frame.ysize: continue

        if not coordinate_box.contains(position): continue

        # Fill the columns
        name_column.append(gal_name)
        #ra_column.append(position.ra.value)
        ra_column.append(position.ra)
        #dec_column.append(position.dec.value)
        dec_column.append(position.dec)
        redshift_column.append(gal_redshift)
        type_column.append(gal_type)
        #alternative_names_column.append(", ".join(gal_names) if len(gal_names) > 0 else None)
        alternative_names_column.append(gal_names)
        #distance_column.append(gal_distance.value if gal_distance is not None else None)
        distance_column.append(gal_distance)
        #inclination_column.append(gal_inclination.degree if gal_inclination is not None else None)
        inclination_column.append(gal_inclination)
        #d25_column.append(gal_d25.value if gal_d25 is not None else None)
        d25_column.append(gal_d25)
        #major_column.append(gal_major.value if gal_major is not None else None)
        major_column.append(gal_major)
        #minor_column.append(gal_minor.value if gal_minor is not None else None)
        minor_column.append(gal_minor)
        #pa_column.append(gal_pa.degree if gal_pa is not None else None)
        pa_column.append(gal_pa)

    # Determine the number of galaxies in the lists
    number_of_galaxies = len(name_column)

    # Indicate which galaxy is the principal galaxy
    principal_column = [False] * number_of_galaxies

    # Determine the index of the galaxy with the maximal major axis length
    principal_index = max(range(number_of_galaxies), key=lambda index: (major_column[index] if major_column[index] is not None else 0.0))
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

    # Return the data
    return name_column, ra_column, dec_column, redshift_column, type_column, alternative_names_column, distance_column, \
           inclination_column, d25_column, major_column, minor_column, pa_column, principal_column, companions_column, \
           parent_column

# -----------------------------------------------------------------

def get_galaxy_s4g_one_component_info(name):

    """
    This function ...
    :param name:
    :return:
    """

    # The Vizier querying object
    vizier = Vizier()
    vizier.ROW_LIMIT = -1

    # Get the "galaxies" table
    result = vizier.query_object(name, catalog=["J/ApJS/219/4/galaxies"])

    # No results?
    if len(result) == 0: return None, None, None, None, None, None

    # Get table
    table = result[0]

    # PA: [0.2/180] Outer isophote position angle
    # e_PA: [0/63] Standard deviation in PA
    # Ell:  [0.008/1] Outer isophote ellipticity
    # e_Ell: [0/0.3] Standard deviation in Ell

    # PA1: Elliptical isophote position angle in deg
    # n: Sersic index
    # Re: effective radius in arcsec

    # Tmag: total magnitude

    s4g_name = table["Name"][0]

    pa = Angle(table["PA"][0] - 90., "deg")
    pa_error = Angle(table["e_PA"][0], "deg")

    ellipticity = table["Ell"][0]
    ellipticity_error = table["e_Ell"][0]

    n = table["n"][0]

    re = table["Re"][0] * u("arcsec")

    mag = table["Tmag"][0]

    # Return the results
    return s4g_name, pa, ellipticity, n, re, mag

# -----------------------------------------------------------------

def get_ra_dec_degrees_from_table(table, index):

    """
    This function ...
    :param table: 
    :param index: 
    :return: 
    """

    # Get units
    ra_unit = table["RAJ2000"].unit
    dec_unit = table["DEJ2000"].unit

    # Get values
    ra_value = table["RAJ2000"][index]
    dec_value = table["DEJ2000"][index]

    # Combine
    return get_ra_dec_degrees_from_values(ra_value, dec_value, ra_unit, dec_unit)

# -----------------------------------------------------------------

def get_ra_dec_degrees_from_values(ra, dec, ra_unit, dec_unit):

    """
    This function ...
    :param ra:
    :param dec:
    :param ra_unit:
    :param dec_unit:
    :return:
    """

    from ...core.tools.strings import unquote
    from . import coordinates

    if isinstance(ra_unit, UnrecognizedUnit):

        ra_unit_string = unquote(str(ra_unit.name))
        if ra_unit_string == "h:m:s":

            ra_degrees = float(coordinates.hms_to_degrees(ra=ra))
            ra_with_unit = ra_degrees * u("deg")

        else: raise ValueError("The only allowed unit in string format for RA is 'h:m:s'")

    else: ra_with_unit = ra * ra_unit

    if isinstance(dec_unit, UnrecognizedUnit):

        dec_unit_string = unquote(str(dec_unit.name))
        if dec_unit_string == "d:m:s":

            dec_degrees = float(coordinates.hms_to_degrees(dec=dec))
            dec_with_unit = dec_degrees * u("deg")

        else: raise ValueError("The only allowed unit in string format for DEC is 'd:m:s'")

    else: dec_with_unit = dec * dec_unit

    #print(type(ra_with_unit), type(dec_with_unit))

    # Return in degrees
    return ra_with_unit.to("deg").value, dec_with_unit.to("deg").value

# -----------------------------------------------------------------

def get_galaxy_info(name, position=None):

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

        #print(ned_entry)

        # Get a more common name for this galaxy (sometimes, the name obtained from NED is one starting with 2MASX .., use the PGC name in this case)
        if ned_entry["Object Name"].startswith("2MASX "): gal_name = name
        else: gal_name = ned_entry["Object Name"]

        # Get the redshift
        gal_redshift = ned_entry["Redshift"]
        if isinstance(gal_redshift, np.ma.core.MaskedConstant): gal_redshift = None

        # Get the type (G=galaxy, HII ...)
        gal_type = ned_entry["Type"]
        if isinstance(gal_type, np.ma.core.MaskedConstant): gal_type = None

        # Get the distance
        # THIS IS PROBABLY NOT AN ACTUAL DISTANCE?
        #ned_distance = ned_entry["Distance (arcmin)"]
        #if isinstance(ned_distance, np.ma.core.MaskedConstant): ned_distance = None
        ned_distance = None

    except astroquery.exceptions.RemoteServiceError:

        # Set attributes
        gal_name = name
        gal_redshift = None
        gal_type = None
        ned_distance = None

    # Create a new Vizier object and set the row limit to -1 (unlimited)
    viz = Vizier(keywords=["galaxies", "optical"])
    viz.ROW_LIMIT = -1

    # Query Vizier and obtain the resulting table
    result = viz.query_object(name.replace(" ", ""), catalog=["VII/237"])

    # Not found ... TODO: fix this ... this object was in the first query output
    if len(result) == 0: return name, position, None, None, [], None, None, None, None, None, None

    table = result[0]

    #print(type(table["RAJ2000"].unit))
    #print(table["_RAJ2000"].unit)
    #print(table["_DEJ2000"].unit)
    #print(table["RAJ2000"].unit)
    #print(table["DEJ2000"].unit)

    ra_unit = table["RAJ2000"].unit
    dec_unit = table["DEJ2000"].unit

    # Get the correct entry (sometimes, for example for mergers, querying with the name of one galaxy gives two hits! We have to obtain the right one each time!)
    if len(table) == 0: raise ValueError("The galaxy could not be found under this name")
    elif len(table) == 1:
        #entry_index = 0
        entry = table[0]
    else:

        #entry_index = None
        entry = None

        # Some rows don't have names, if no match is found based on the name just take the row that has other names defined
        rows_with_names = []
        for row in table:
            if row["ANames"]: rows_with_names.append(row)

        # If only one row remains, take that one for the galaxy we are looking for
        if len(rows_with_names) == 1: entry = rows_with_names[0]

        # Else, loop over the rows where names are defined and look for a match
        else:
            for row in rows_with_names:

                names = row["ANames"]

                if name.replace(" ", "") in names or gal_name.replace(" ", "") in names:

                    entry = row
                    break

        # If no matches are found, look for the table entry for which the coordinate matches the given position (if any)
        if entry is None and position is not None:
            #for row in table:
            for index in range(len(table)):
                ra, dec = get_ra_dec_degrees_from_table(table, index)
                if np.isclose(ra, position.ra.to("deg").value) and np.isclose(dec, position.dec.to("deg").value):
                    entry = row
                    break

    # Note: another temporary fix
    if entry is None: return name, position, None, None, [], None, None, None, None, None, None

    # Get coordinate
    ra, dec = get_ra_dec_degrees_from_values(entry["RAJ2000"], entry["DEJ2000"], ra_unit, dec_unit)

    # Get the right ascension and the declination
    #position = SkyCoordinate(ra=entry["RAJ2000"], dec=entry["DEJ2000"], unit="deg", frame="fk5")
    position = SkyCoordinate(ra=ra, dec=dec, unit="deg", frame="fk5")

    # Get the names given to this galaxy
    gal_names = entry["ANames"].split() if entry["ANames"] else []

    # Get the size of the galaxy
    ratio = np.power(10.0, entry["logR25"]) if entry["logR25"] else None
    diameter = np.power(10.0, entry["logD25"]) * 0.1 * u("arcmin") if entry["logD25"] else None

    #print("  ratio = ", ratio)
    #print("  D25_diameter = ", diameter)
    #print("  position = ", position)

    radial_profiles_result = viz.query_object(name, catalog="J/ApJ/658/1006")

    if len(radial_profiles_result) > 0:

        radial_profiles_entry = radial_profiles_result[0][0]

        gal_distance = radial_profiles_entry["Dist"] * u("Mpc")
        gal_inclination = Angle(radial_profiles_entry["i"], "deg")
        gal_d25 = radial_profiles_entry["D25"] * u("arcmin")

    else:

        gal_distance = ned_distance
        gal_inclination = None
        gal_d25 = diameter

    # Distance not found?
    if gal_distance is None:

        grav_result = viz.query_object(name, catalog="VII/267/gwgc")
        if len(grav_result) > 0:

            table = grav_result[0]
            gal_distance = table["Dist"][0] * u("Mpc")

    # Get the size of major and minor axes
    gal_major = diameter
    gal_minor = diameter / ratio if diameter is not None and ratio is not None else None

    #print(" gal_major", diameter)
    #print(" gal_minor", gal_minor)

    # Get the position angle of the galaxy
    gal_pa = Angle(entry["PA"] - 90.0, "deg") if entry["PA"] else None

    #print(" gal_pa", gal_pa)

    # Create and return a new Galaxy instance
    return gal_name, position, gal_redshift, gal_type, gal_names, gal_distance, gal_inclination, gal_d25, gal_major, gal_minor, gal_pa

# -----------------------------------------------------------------

def galaxies_in_box(center, ra_span, dec_span):

    """
    This function ...
    :param center:
    :param ra_span:
    :param dec_span:
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

    # Debugging
    log.debug("Querying the HYPERLEDA catalog ...")

    # Query Vizier and obtain the resulting table
    result = viz.query_region(center.to_astropy(), width=ra_span, height=dec_span, catalog=["VII/237"])

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
    if result is None or len(result) == 0:

        ra_span *= 1.0 + 1e-5
        result = viz.query_region(center.to_astropy(), width=ra_span, height=dec_span, catalog=["VII/237"])

    table = result[0]

    # Loop over the rows in the table
    for entry in table:
        name = "PGC " + str(entry["PGC"])
        coordinate = SkyCoordinate(ra=entry["RAJ2000"], dec=entry["DEJ2000"], unit="deg", frame="fk5")
        namepluscoordinate = (name, coordinate)
        names.append(namepluscoordinate)

    # Return the list of galaxies
    return names

# -----------------------------------------------------------------
