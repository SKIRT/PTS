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
from ...core.basics.log import log
from ...core.tools import filesystem as fs
from ...core.launch.pts import PTSRemoteLauncher
from ...core.tools import network, archive
from ...magic.core.frame import Frame
from ...core.tools.serialization import write_dict
from ...core.launch.pts import launch_local
from .component import galex, sdss, twomass, spitzer, wise, herschel, planck, other, halpha

# -----------------------------------------------------------------

class ImageFetcher(DataComponent):

    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(ImageFetcher, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The DustPedia database
        self.database = DustPediaDatabase()

        # The urls of the images found on the DustPedia archive, for each origin
        self.dustpedia_image_urls = defaultdict(dict)

        # Create the PTS remote environment
        self.launcher = None

    # -----------------------------------------------------------------

    @property
    def has_halpha_url(self):
        return self.config.halpha_url is not None

    # -----------------------------------------------------------------

    @property
    def has_other_urls(self):
        return self.config.other_urls is not None and len(self.config.other_urls) > 0

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Fetch the images urls from the DustPedia archive
        self.get_dustpedia_urls()

        # 3. Fetch GALEX data and calculate poisson errors
        self.fetch_galex()

        # 4. Fetch SDSS data and calculate poisson errors
        self.fetch_sdss()

        # 5. Fetch the H-alpha image
        if self.has_halpha_url: self.fetch_halpha()

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

        # Fetch other images
        if self.has_other_urls: self.fetch_other()

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

            # Loop over all URLs, indexed on image name
            for name in all_urls:

                # Skip error frames unless the 'errors' flag has been enabled
                if not self.config.errors and "_Error" in name: continue

                # Add url to the dictionary
                if origin == herschel:
                    if "pacs" in name.lower() or "spire" in name.lower(): self.dustpedia_image_urls[origin][name] = all_urls[name]

                # Not Herschel
                elif origin in name: self.dustpedia_image_urls[origin][name] = all_urls[name]

    # -----------------------------------------------------------------

    def fetch_from_dustpedia(self, origin, common_origin=None):

        """
        This function ...
        :return:
        """

        if common_origin is None: common_origin = origin

        # Loop over all images from this origin
        for name in self.dustpedia_image_urls[origin]:

            # Determine the path to the image file
            path = fs.join(self.data_images_paths[common_origin], name)

            # Check if the image is already present
            if fs.is_file(path):
                log.success("The '" + name + "' image is already present")
                continue

            # Debugging
            log.debug("Fetching the '" + name + "' image from the DustPedia archive ...")

            # Download the image
            url = self.dustpedia_image_urls[origin][name]
            self.database.download_image_from_url(url, path)

            # Success
            log.success("The '" + name + "' image was downloaded")

    # -----------------------------------------------------------------

    def fetch_galex(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching the GALEX images ...")

        # Fetch the GALEX data from the DustPedia archive
        self.fetch_from_dustpedia(galex)

        # Make the GALEX poisson error maps
        if self.config.make_poisson: self.make_poisson_galex()
        else: log.warning("The GALEX poisson error maps will have to be created manually with 'make_galex' and be placed next to the images with the suffix '_poisson'")

    # -----------------------------------------------------------------

    def make_poisson_galex(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Launching the procedures to create GALEX poisson error mosaic maps ...")

        # Determine local output directory path
        local_output_path = fs.create_directory_in(self.data_images_paths[galex], "temp")

        # Create the configuration dictionary
        config_dict = dict()
        config_dict["galaxy_name"] = self.ngc_name_nospaces
        config_dict["output"] = local_output_path
        config_dict["max_nobservations_fuv"] = self.config.max_nobservations_mosaic
        config_dict["max_nobservations_nuv"] = self.config.max_nobservations_mosaic
        config_dict["nprocesses"] = self.config.nprocesses

        # Set the analysis info and analyser class
        analysis_info = {"modeling_path": self.config.path}
        analysers = ["pts.modeling.data.analyser.MosaicAnalyser"]

        # Create the GALEX mosaic and Poisson errors frame
        command = "make_galex"

        # Local
        if self.launcher is None: launch_local(command, config_dict, analysers=analysers, analysis_info=analysis_info)

        # Remote, attached
        elif self.config.attached: self.launcher.run_and_analyse(command, config_dict, local_output_path, analysers, analysis_info, use_session=False)

        # Remote, run in detached mode
        else:
            self.launcher.run_detached("make_galex", config_dict, analysers=analysers, analysis_info=analysis_info, remove_local_output=True)
            self.detached = True

    # -----------------------------------------------------------------

    def fetch_sdss(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching the SDSS images ...")

        # Fetch the SDSS data from the DustPedia archive
        self.fetch_from_dustpedia(sdss)

        # Make the SDSS poisson error maps
        if self.config.make_poisson: self.make_poisson_sdss()
        else: log.warning("The SDSS poisson error maps will have to be created manually with 'make_sdss' and be placed next to the images with the suffix '_poisson'")

    # -----------------------------------------------------------------

    def make_poisson_sdss(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Launching the procedures to create SDSS poisson error mosaic maps ...")

        # Determine the output directory path
        local_output_path = fs.create_directory_in(self.data_images_paths[sdss], "temp")

        # Create the configuration dictionary
        config_dict = dict()
        config_dict["galaxy_name"] = self.ngc_name_nospaces
        config_dict["output"] = fs.join(local_output_path)
        config_dict["max_nobservations_u"] = self.config.max_nobservations_mosaic
        config_dict["max_nobservations_g"] = self.config.max_nobservations_mosaic
        config_dict["max_nobservations_r"] = self.config.max_nobservations_mosaic
        config_dict["max_nobservations_i"] = self.config.max_nobservations_mosaic
        config_dict["max_nobservations_z"] = self.config.max_nobservations_mosaic
        config_dict["nprocesses"] = self.config.nprocesses

        # Set the analysis info and analyser class
        analysis_info = {"modeling_path": self.config.path}
        analysers = ["pts.modeling.data.analyser.MosaicAnalyser"]

        # Create the SDSS mosaic and Poisson errors frame
        command = "make_sdss"

        # Local
        if self.launcher is None: launch_local(command, config_dict, analysers=analysers, analysis_info=analysis_info)

        # Remote attached
        elif self.config.attached: self.launcher.run_and_analyse(command, config_dict, local_output_path, analysers, analysis_info, use_session=False)

        # Run in detached mode
        else:
            self.launcher.run_detached(command, config_dict, analysers=analysers, analysis_info=analysis_info, remove_local_output=True)
            self.detached = True

    # -----------------------------------------------------------------

    def fetch_halpha(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching the H-alpha image ...")

        # Download the Halpha image
        image_path = network.download_file(self.config.halpha_url, self.data_images_paths[halpha])

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
        self.fetch_from_dustpedia(twomass)

    # -----------------------------------------------------------------

    def fetch_spitzer(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching the Spitzer images ...")

        # Fetch the Spitzer data from the DustPedia archive
        self.fetch_from_dustpedia(spitzer)

    # -----------------------------------------------------------------

    def fetch_wise(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching the WISE images ...")

        # Fetch the WISE data from the DustPedia archive
        self.fetch_from_dustpedia(wise)

    # -----------------------------------------------------------------

    def fetch_herschel(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching the Herschel images ...")

        # Fetch the Pacs data
        #self.fetch_from_dustpedia("PACS", "Herschel")

        # Fetch SPIRE
        #self.fetch_from_dustpedia("SPIRE", "Herschel")

        self.fetch_from_dustpedia(herschel)

    # -----------------------------------------------------------------

    def fetch_planck(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching the Planck images ...")

        # Fetch the Planck data from the DustPedia archive
        self.fetch_from_dustpedia(planck)

    # -----------------------------------------------------------------

    def fetch_other(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Fetching other images ...")

        # Download the images
        paths = network.download_files(self.config.other_urls, self.data_images_path[other])

        # Unpack the image files if necessary
        for path in paths:
            if not path.endswith("fits"): path = archive.decompress_file_in_place(path, remove=True)

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

        # Write
        write_dict(self.dustpedia_image_urls, self.urls_path)

# -----------------------------------------------------------------
