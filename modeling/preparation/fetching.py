#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.preparation.fetching Contains the DataFetcher class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import PreparationComponent
from ...magic.misc.dustpedia import DustPediaDatabase
from ...core.tools.logging import log

# -----------------------------------------------------------------

class DataFetcher(PreparationComponent):
    
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

        # The DustPedia database
        self.database = DustPediaDatabase()

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

    # -----------------------------------------------------------------

    def fetch_galaxy_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching galaxy info ...")

        # Get the info
        info = self.database.get_galaxy_info(self.galaxy_name)

        print(info)

    # -----------------------------------------------------------------

    def fetch_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching the images ...")

        # Get the image names
        image_names = self.database.get_image_names(self.galaxy_name)

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
        sed = self.database.get_sed(self.galaxy_name)

# -----------------------------------------------------------------
