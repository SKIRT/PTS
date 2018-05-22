#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.dustpedia.data.images Contains the ImageFetcher class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.basics.log import log
from ...core.tools import filesystem as fs
from ..core.sample import DustPediaSample
from ..core.database import DustPediaDatabase, get_account

# -----------------------------------------------------------------

instruments = ['2MASS', 'DSS', 'GALEX', 'PACS', 'Planck', 'SDSS', 'SPIRE', 'Spitzer', 'WISE']

# -----------------------------------------------------------------

class ImageFetcher(Configurable):
        
    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ImageFetcher, self).__init__(*args, **kwargs)

        # The DustPediaSample object
        self.sample = DustPediaSample()

        # The DustPediaDataBase object
        self.database = DustPediaDatabase()

        # The NGC name
        self.ngc_name = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Get the images
        self.get_images()

        # Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ImageFetcher, self).setup(**kwargs)

        # Get the NGC name
        self.ngc_name = self.sample.get_name(self.config.galaxy_name)

        # Get username and password for the DustPedia database
        if self.config.database.username is not None:
            username = self.config.database.username
            password = self.config.database.password
        else: username, password = get_account()

        # Login to the DustPedia database
        self.database.login(username, password)

    # -----------------------------------------------------------------

    def get_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the images ...")

        #filters = self.database.get_image_filters(self.ngc_name)
        names_and_urls = self.database.get_image_names_and_urls(self.ngc_name)

        paths = dict()

        # Create directories
        for instrument in self.config.instruments:

            # Inform the user
            log.info("Creating directory for the " + instrument + " instrument ...")

            # Create directory
            path = fs.create_directory_in(self.config.path, instrument)
            paths[instrument] = path

        # Loop over all entries
        for name, url in names_and_urls.items():

            # Get the instrument
            instrument = url.split("instrument=")[1]

            # Skip if not specified
            if instrument not in self.config.instruments: continue

            # Determine image path
            path = fs.join(paths[instrument], name)

            # Inform the user
            log.info("Downloading the " + name + " image ...")

            # Download via database (decompressing is also done)
            self.database.download_image_from_url(url, path)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------
