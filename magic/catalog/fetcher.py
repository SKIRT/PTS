#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.catalog.fetcher Contains the CatalogFetcher class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..tools import catalogs
from ...core.basics.log import log
from .extended import ExtendedSourceCatalog
from .point import PointSourceCatalog
from ...core.tools import types
from ...core.tools.stringify import tostr

# -----------------------------------------------------------------

class CatalogFetcher(object):

    """
    This function ...
    """

    def __init__(self):

        """
        The constructor ...
        """

        pass

    # -----------------------------------------------------------------

    @staticmethod
    def get_extended_source_catalog(coordinate_box):

        """
        This function ...
        :param coordinate_box:
        :return:
        """

        # Inform the user
        log.info("Searching extended sources in: ")
        log.info(" - center: " + tostr(coordinate_box.center))
        log.info(" - radius: " + tostr(coordinate_box.radius))

        # Get the data
        name_column, ra_column, dec_column, redshift_column, type_column, alternative_names_column, distance_column, \
        inclination_column, d25_column, major_column, minor_column, pa_column, principal_column, companions_column, \
        parent_column = catalogs.create_galaxy_catalog(coordinate_box)

        # Create catalog
        catalog = ExtendedSourceCatalog()

        # Add entries
        for index in range(len(name_column)):

            catalog.add_entry(name_column[index], ra_column[index], dec_column[index], redshift_column[index], type_column[index],
                              alternative_names_column[index], distance_column[index], inclination_column[index],
                              d25_column[index], major_column[index], minor_column[index], pa_column[index],
                              principal_column[index], companions_column[index], parent_column[index])

        # Return the catalog
        return catalog

    # -----------------------------------------------------------------

    @staticmethod
    def get_point_source_catalog(coordinate_box, pixelscale, catalogues):

        """
        This function ...
        :param coordinate_box:
        :param pixelscale:
        :param catalogues:
        :return:
        """

        # Inform the user
        log.info("Searching point sources in: ")
        log.info(" - center: " + tostr(coordinate_box.center))
        log.info(" - radius: " + tostr(coordinate_box.radius))

        # Check whether the 'catalogs' setting defines a single catalog name or a list of such names
        if types.is_string_type(catalogues): catalog_list = [catalogues]
        elif isinstance(catalogues, list): catalog_list = catalogues
        else: raise ValueError("Invalid option for 'catalogs', should be a string or a list of strings")

        # Get the data
        catalog_column, id_column, ra_column, dec_column, ra_error_column, dec_error_column, \
        confidence_level_column = catalogs.create_star_catalog(coordinate_box, pixelscale, catalog_list)

        # Create catalog
        catalog = PointSourceCatalog()

        # Add entries
        for index in range(len(catalog_column)):

            # Debugging
            log.info("Adding entry " + str(index+1) + " to the catalog ...")

            catalog.add_entry(catalog_column[index], id_column[index], ra_column[index], dec_column[index],
                              ra_error_column[index], dec_error_column[index], confidence_level_column[index])

        # Return the catalog
        return catalog

# -----------------------------------------------------------------
