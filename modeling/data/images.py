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

# Import standard modules
from collections import defaultdict

# Import the relevant PTS classes and modules
from .component import DataComponent
from ...dustpedia.core.database import DustPediaDatabase, get_account
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
from ...core.launch.pts import PTSRemoteLauncher

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

        # The names of the images found on the DustPedia archive, for each observatory
        self.dustpedia_image_names = defaultdict(list)

        # Create the PTS remote environment
        self.launcher = PTSRemoteLauncher()

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Fetch the images from the DustPedia archive
        self.get_dustpedia_names()

        # 3. Fetch GALEX data and calculate poisson errors
        self.fetch_galex()

        # 4. Fetch SDSS data and calculate poisson errors
        self.fetch_sdss()

        # 5. Fetch the H-alpha image
        #self.fetch_halpha()

        # 6. Fetch the 2MASS images
        #self.fetch_2mass()

        # 7. Fetch the Spitzer images
        #self.fetch_spitzer()

        # 8. Fetch the WISE images
        #self.fetch_wise()

        # 9. Fetch the Herschel images
        #self.fetch_herschel()

        # 10. Writing
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

        # Setup the remote PTS launcher
        self.launcher.setup(self.config.remote)

    # -----------------------------------------------------------------

    def get_dustpedia_names(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching the names of the images on the DustPedia database ...")

        # Get the image names
        all_names = self.database.get_image_names(self.ngc_id_nospaces)

        # Order the names per origin
        for origin in self.data_origins:
            for name in all_names:

                if not self.config.errors and "_Error" in name: continue # Skip error frames unless the 'errors' flag has been enabled
                if origin in name: self.dustpedia_image_names[origin].append(name)

    # -----------------------------------------------------------------

    def fetch_from_dustpedia(self, origin):

        """
        This function ...
        :return:
        """

        # Loop over all images from this origin
        for name in self.dustpedia_image_names[origin]:

            # Determine the path to the image file
            path = fs.join(self.data_images_paths[origin], name)

            # Download the image
            self.database.download_image(self.ngc_id_nospaces, name, path)

    # -----------------------------------------------------------------

    def fetch_galex(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching the GALEX data ...")

        # Fetch the GALEX data from the DustPedia archive
        self.fetch_from_dustpedia("GALEX")

        # Create the configuration dictionary
        config_dict = dict()
        config_dict["galaxy_name"] = self.ngc_id_nospaces
        config_dict["output"] = fs.join(self.data_images_paths["GALEX"], "temp")
        fs.create_directory(config_dict["output"])

        # Set the analysis info and analyser class
        analysis_info = {"modeling_path": self.config.path}
        analysers = ["pts.modeling.data.analyser.MosaicAnalyser"]

        # Create the GALEX mosaic and Poisson errors frame
        self.launcher.run_detached("make_galex", config_dict, analysers=analysers, analysis_info=analysis_info)

    # -----------------------------------------------------------------

    def fetch_sdss(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching the SDSS data ...")

        # Fetch the SDSS data from the DustPedia archive
        self.fetch_from_dustpedia("SDSS")

        # Create the configuration dictionary
        config_dict = dict()
        config_dict["galaxy_name"] = self.ngc_id_nospaces
        config_dict["output"] = fs.join(self.data_images_paths["SDSS"], "temp")
        fs.create_directory(config_dict["output"])

        # Set the analysis info and analyser class
        analysis_info = {"modeling_path": self.config.path}
        analysers = ["pts.modeling.data.analyser.MosaicAnalyser"]

        # Create the SDSS mosaic and Poisson errors frame
        self.launcher.run_detached("make_sdss", config_dict, analysers=analysers, analysis_info=analysis_info)

    # -----------------------------------------------------------------

    def fetch_halpha(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching the H-alpha image ...")

    # -----------------------------------------------------------------

    def fetch_2mass(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching the ...")

        # Fetch the 2MASS data from the DustPedia archive
        self.fetch_from_dustpedia("2MASS")

    # -----------------------------------------------------------------

    def fetch_spitzer(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching the Spitzer data ...")

        # Fetch the Spitzer data from the DustPedia archive
        self.fetch_from_dustpedia("Spitzer")

    # -----------------------------------------------------------------

    def fetch_wise(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching the ...")

        # Fetch the WISE data from the DustPedia archive
        self.fetch_from_dustpedia("WISE")

    # -----------------------------------------------------------------

    def fetch_herschel(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching the ...")

        # Fetch the Herschel data from the DustPedia archive
        self.fetch_from_dustpedia("Herschel")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

# -----------------------------------------------------------------
