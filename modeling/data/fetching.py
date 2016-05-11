#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.data.fetching Contains the DataFetcher class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import DataComponent
from ...magic.misc.dustpedia import DustPediaDatabase
from ...core.tools.logging import log
from ...magic.tools import catalogs
from ...core.tools import tables

# -----------------------------------------------------------------

class DataFetcher(DataComponent):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(DataFetcher, self).__init__(config)

        # -- Attributes --

        # The NGC id of the galaxy
        self.ngc_id = None

        # The DustPedia database
        self.database = DustPediaDatabase()

        # The galaxy info
        self.info = None

        # The images
        self.images = []

        # The SED
        self.sed = None

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new DataFetcher instance
        fetcher = cls()

        # Set the modeling path
        fetcher.config.path = arguments.path

        # Return the data fetcher
        return fetcher

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # Get the galaxy info
        self.fetch_galaxy_info()

        # 2. Fetch the images
        self.fetch_images()

        # 3. Fetch the SED
        self.fetch_sed()

        # 4. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(DataFetcher, self).setup()

        # Login to the DustPedia database
        self.database.login(self.config.database.username, self.config.database.password)

        # Determine the NGC id of the galaxy
        self.ngc_id = catalogs.get_ngc_name(self.galaxy_name, delimiter="")

    # -----------------------------------------------------------------

    def fetch_galaxy_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching galaxy info ...")

        # Get the info
        self.info = self.database.get_galaxy_info(self.ngc_id)

    # -----------------------------------------------------------------

    def fetch_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching the images ...")

        # Get the image names
        image_names = self.database.get_image_names(self.ngc_id)

        print(image_names)

    # -----------------------------------------------------------------

    def fetch_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching the SED ...")

        # Get the SED
        self.sed = self.database.get_sed(self.ngc_id)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the galaxy info
        self.write_info()

        # Write the images
        self.write_images()

        # Write the SED
        self.write_sed()

    # -----------------------------------------------------------------

    def write_info(self):

        """
        This function ...
        :return:
        """

        # Infom the user
        log.info("Writing the galaxy info ...")

        # Write the galaxy info table
        tables.write(self.info, self.galaxy_info_path, format="ascii.ecsv")

    # -----------------------------------------------------------------

    def write_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the galaxy images ...")

    # -----------------------------------------------------------------

    def write_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the SED ...")

# -----------------------------------------------------------------
