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
from ..basics.catalogcoverage import CatalogCoverage
from ..tools import catalogs
from ...core.basics.configurable import Configurable
from ...core.tools import introspection, tables
from ...core.tools import filesystem as fs
from ...core.tools.logging import log

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

        catalog = catalogs.create_galaxy_catalog(coordinate_box)

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

        # Check whether the 'catalogs' setting defines a single catalog name or a list of such names
        if isinstance(catalogues, basestring): catalog_list = [catalogues]
        elif isinstance(catalogues, list): catalog_list = catalogues
        else: raise ValueError("Invalid option for 'catalogs', should be a string or a list of strings")

        # Create the star catalog
        catalog = catalogs.create_star_catalog(coordinate_box, pixelscale, catalog_list)

        # Return the catalog
        return catalog

# -----------------------------------------------------------------
