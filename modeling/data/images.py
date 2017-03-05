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
from ...core.tools import network, archive
from ...magic.core.frame import Frame
from ...core.tools.serialization import write_dict
from ...core.launch.pts import launch_local

# -----------------------------------------------------------------

class ImageFetcher(DataComponent):
    
    """
    This class...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(ImageFetcher, self).__init__(config, interactive)

        # -- Attributes --

        # The DustPedia database
        self.database = DustPediaDatabase()

        # The urls of the images found on the DustPedia archive, for each origin
        self.dustpedia_image_urls = defaultdict(dict)

        # Create the PTS remote environment
        self.launcher = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Fetch the images urls from the DustPedia archive
        self.get_dustpedia_urls()

        # 3. Fetch GALEX data and calculate poisson errors
        self.fetch_galex()

        # 4. Fetch SDSS data and calculate poisson errors
        self.fetch_sdss()

        # 5. Fetch the H-alpha image
        self.fetch_halpha()

        # 6. Fetch the 2MASS images
        self.fetch_2mass()

        # 7. Fetch the Spitzer images
        self.fetch_spitzer()

        # 8. Fetch the WISE images
        self.fetch_wise()

        # 9. Fetch the Herschel images
        self.fetch_herschel()

        # 10. Fetch the Planck images
        self.fetch_planck()

        # 11. Writing
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

        # Get username and password for the DustPedia database
        if self.config.database.username is not None:
            username = self.config.database.username
            password = self.config.database.password
        else: username, password = get_account()

        # Login to the DustPedia database
        self.database.login(username, password)

        # Setup the remote PTS launcher
        if self.config.remote is not None:
            self.launcher = PTSRemoteLauncher()
            self.launcher.setup(self.config.remote)

    # -----------------------------------------------------------------

    def get_dustpedia_urls(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching the names of the images on the DustPedia database ...")

        # Get the image names
        all_urls = self.database.get_image_names_and_urls(self.ngc_name_nospaces)

        # Order the names per origin
        for origin in self.data_origins:

            # Loop over all URLs
            for name in all_urls:

                if not self.config.errors and "_Error" in name: continue # Skip error frames unless the 'errors' flag has been enabled
                if origin in name: self.dustpedia_image_urls[origin][name] = all_urls[name]

    # -----------------------------------------------------------------

    def fetch_from_dustpedia(self, origin):

        """
        This function ...
        :return:
        """

        # Loop over all images from this origin
        for name in self.dustpedia_image_urls[origin]:

            # Debugging
            log.debug("Fetching the '" + name + "' image from the DustPedia archive ...")

            # Determine the path to the image file
            path = fs.join(self.data_images_paths[origin], name)

            # Check if the image is already present
            if fs.is_file(path):
                log.warning("The '" + name + "' image is already present")
                continue

            # Download the image
            url = self.dustpedia_image_urls[origin][name]
            self.database.download_image_from_url(url, path)

    # -----------------------------------------------------------------

    def fetch_galex(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching the GALEX images ...")

        # Fetch the GALEX data from the DustPedia archive
        self.fetch_from_dustpedia("GALEX")

        # Determine local output directory path
        local_output_path = fs.create_directory_in(self.data_images_paths["GALEX"], "temp")

        # Create the configuration dictionary
        config_dict = dict()
        config_dict["galaxy_name"] = self.ngc_name_nospaces
        config_dict["output"] = local_output_path

        # Set the analysis info and analyser class
        analysis_info = {"modeling_path": self.config.path}
        analysers = ["pts.modeling.data.analyser.MosaicAnalyser"]

        # Create the GALEX mosaic and Poisson errors frame
        command = "make_galex"

        # Local
        if self.launcher is None: launch_local(command, config_dict, analysers=analysers, analysis_info=analysis_info)

        # Remote, attached
        elif self.config.attached: self.launcher.run_and_analyse(command, config_dict, local_output_path, analysers, analysis_info)

        # Remote, run in detached mode
        else: self.launcher.run_detached("make_galex", config_dict, analysers=analysers, analysis_info=analysis_info, remove_local_output=True)

    # -----------------------------------------------------------------

    def fetch_sdss(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching the SDSS images ...")

        # Fetch the SDSS data from the DustPedia archive
        self.fetch_from_dustpedia("SDSS")

        # Determine the output directory path
        local_output_path = fs.create_directory_in(self.data_images_paths["SDSS"], "temp")

        # Create the configuration dictionary
        config_dict = dict()
        config_dict["galaxy_name"] = self.ngc_name_nospaces
        config_dict["output"] = fs.join(local_output_path)

        # Set the analysis info and analyser class
        analysis_info = {"modeling_path": self.config.path}
        analysers = ["pts.modeling.data.analyser.MosaicAnalyser"]

        # Create the SDSS mosaic and Poisson errors frame
        command = "make_sdss"

        # Local
        if self.launcher is None: launch_local(command, config_dict, analysers=analysers, analysis_info=analysis_info)

        # Remote attached
        elif self.config.attached: self.launcher.run_and_analyse(command, config_dict, local_output_path, analysers, analysis_info)

        # Run in detached mode
        else: self.launcher.run_detached(command, config_dict, analysers=analysers, analysis_info=analysis_info, remove_local_output=True)

    # -----------------------------------------------------------------

    def fetch_halpha(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching the H-alpha image ...")

        # Download the Halpha image
        image_path = network.download_file(self.config.halpha_url, self.data_images_paths["Halpha"])

        # Unpack the image file if necessary
        if not image_path.endswith("fits"): image_path = archive.decompress_file_in_place(image_path, remove=True)

        # Rescale the image to have a certain flux value, if specified
        if self.config.halpha_flux is not None:

            # Inform the user
            log.info("Rescaling the H-alpha image to a flux value of " + str(self.config.halpha_flux) + " ...")

            # Open the image
            frame = Frame.from_file(image_path)

            # Normalize to the flux value
            frame.normalize(self.config.halpha_flux.value)

            # Set the unit
            frame.unit = self.config.halpha_flux.unit

            # Save the image
            frame.saveto(image_path)

    # -----------------------------------------------------------------

    def fetch_2mass(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching the 2MASS images ...")

        # Fetch the 2MASS data from the DustPedia archive
        self.fetch_from_dustpedia("2MASS")

    # -----------------------------------------------------------------

    def fetch_spitzer(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching the Spitzer images ...")

        # Fetch the Spitzer data from the DustPedia archive
        self.fetch_from_dustpedia("Spitzer")

    # -----------------------------------------------------------------

    def fetch_wise(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching the WISE images ...")

        # Fetch the WISE data from the DustPedia archive
        self.fetch_from_dustpedia("WISE")

    # -----------------------------------------------------------------

    def fetch_herschel(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching the Herschel images ...")

        # Fetch the Herschel data from the DustPedia archive
        self.fetch_from_dustpedia("Herschel")

    # -----------------------------------------------------------------

    def fetch_planck(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching the Planck images ...")

        # Fetch the Planck data from the DustPedia archive
        self.fetch_from_dustpedia("Planck")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the URLs
        self.write_urls()

    # -----------------------------------------------------------------

    def write_urls(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the image URLs ...")

        # Determine the path
        path = fs.join(self.data_images_path, "urls.dat")

        # Write
        write_dict(self.dustpedia_image_urls, path)

# -----------------------------------------------------------------
