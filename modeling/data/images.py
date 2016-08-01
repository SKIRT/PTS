#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.data.images Contains the ImageFetcher class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import DataComponent
from ...dustpedia.core.database import DustPediaDatabase, get_account
from ...core.tools.logging import log
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

class ImageFetcher(DataComponent):
    
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
        super(ImageFetcher, self).__init__(config)

        # -- Attributes --

        # The DustPedia database
        self.database = DustPediaDatabase()

        # The images
        self.images = []

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Fetch the images
        self.fetch_images()

        # 3. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(ImageFetcher, self).setup()

        # Get username and password for the DustPedia database
        if self.config.database.username is not None:
            username = self.config.database.username
            password = self.config.database.password
        else: username, password = get_account()

        # Login to the DustPedia database
        self.database.login(username, password)

    # -----------------------------------------------------------------

    def fetch_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching the images ...")

        # Get the image names
        image_names = self.database.get_image_names(self.ngc_id_nospaces)

        # Download the images
        for image_name in image_names:

            # Determine the path to the image file
            path = fs.join(self.data_path, image_name)

            # Download the image
            self.database.download_image(self.ngc_id_nospaces, image_name, path)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

# -----------------------------------------------------------------
