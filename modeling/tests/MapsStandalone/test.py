#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.core.tools.logging import log
from pts.magic.maps.colour.colour import make_map as make_colour_map
from pts.magic.maps.ssfr.colours import make_map as make_ssfr_map
from pts.magic.maps.youngstars.young import make_maps as make_young_stellar_maps
from pts.magic.maps.oldstars.disk import make_map as make_old_stellar_map
from pts.magic.maps.attenuation.cortese import make_map as make_fuv_attenuation_map
from pts.magic.maps.ionizingstars.ionizing import make_map as make_ionizing_stellar_map
from pts.magic.maps.dust.attenuation import make_map
from pts.core.test.implementation import TestImplementation
from pts.dustpedia.core.database import DustPediaDatabase, get_account
from pts.dustpedia.core.sample import resolve_name, get_distance
from pts.core.tools import introspection
from pts.magic.core.list import FrameList
from pts.core.tools import types
from pts.core.filter.filter import parse_filter

# -----------------------------------------------------------------

description = "testing the map making for a certain galaxy"

# -----------------------------------------------------------------

maps_commands = ["make_colours_maps", "make_ssfr_maps", "make_tir_maps", "make_attenuation_maps", "make_old_stars_map", "make_dust_map", "make_young_stars_map", "make_ionizing_stars_map"]

# -----------------------------------------------------------------

class MapsStandaloneTest(TestImplementation):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor
        :param kwargs:
        """

        # Call the constructor of the base class
        super(MapsStandaloneTest, self).__init__(*args, **kwargs)

        # The galaxy distance
        self.distance = None

        # The data path
        self.data_path = None

        # The DustPedia database
        self.database = None

        # The frames
        self.frames = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Setup
        self.setup(**kwargs)

        # 2. Get the data
        if not self.from_existing_reference: self.get_data()
        else: self.check_reference_data()

        # Load the frames
        self.load_frames()

        # Make colour maps
        self.make_colour_maps()

        # Make sSFR maps
        self.make_ssfr_maps()

        # Make TIR maps
        self.make_tir_maps()

        # Make attenuation maps
        self.make_attenuation_maps()

        # Make maps of old stars
        self.make_old_stars_maps()

        # Make dust maps
        self.make_dust_maps()

        # Make young stars maps
        self.make_young_stars_maps()

        # Make ionizing stars maps
        self.make_ionizing_stars_maps()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(MapsStandaloneTest, self).setup(**kwargs)

        # Create data path
        self.data_path = fs.create_directory_in(self.path, "data")

        # Login to the DustPedia database
        self.database = DustPediaDatabase()
        username, password = get_account()
        self.database.login(username, password)

        # Get the distance
        self.distance = get_distance(self.config.galaxy)

    # -----------------------------------------------------------------

    @property
    def from_existing_reference(self):

        """
        This function ...
        :return:
        """

        return self.config.reference_path is not None or self.config.reference_test is not None

    # -----------------------------------------------------------------

    def get_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the image data ...")

        # Resolve the name
        galaxy_name = resolve_name(self.config.galaxy)

        # Download the images
        self.database.download_images(galaxy_name, self.data_path, error_maps=False, not_instruments=["DSS"])

    # -----------------------------------------------------------------

    def check_reference_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the reference data ...")

        # Determine simulation directory
        if self.config.reference_path is not None: data_path = self.config.reference_path
        elif self.config.reference_test is not None: data_path = fs.join(introspection.pts_tests_dir, self.config.reference_test, "data")
        else: raise ValueError("Reference path and reference test settings are None")

        # Check whether directory exist and not empty
        if not fs.is_directory(data_path): raise ValueError("Directory does not exist: " + data_path)
        if fs.is_empty(data_path): raise ValueError("Empty directory: " + data_path)

        # Remove data directory for this test
        fs.remove_directory(self.data_path)

        # Set data path
        self.data_path = data_path

    # -----------------------------------------------------------------

    def load_frames(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the frames ...")

        # Load as frame list
        self.frames = FrameList.from_directory(self.data_path, recursive=True, not_contains="_DSS")

        # Set the distance (to each frame)
        self.frames.distance = self.distance

    # -----------------------------------------------------------------

    def get_frame(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        if types.is_string_type(fltr): fltr = parse_filter(fltr)
        return self.frames[fltr]

    # -----------------------------------------------------------------

    def make_colour_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making colour maps ...")

        # Make FUV-H
        fuv_h = make_colour_map(self.get_frame("FUV"), self.get_frame("H"))

        print(fuv_h)

    # -----------------------------------------------------------------

    def make_ssfr_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making sSFR maps ...")

    # -----------------------------------------------------------------

    def make_tir_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making TIR maps ...")

    # -----------------------------------------------------------------

    def make_attenuation_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making attenuation maps ...")

    # -----------------------------------------------------------------

    def make_old_stars_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making old stellar maps ...")

    # -----------------------------------------------------------------

    def make_dust_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making dust maps ...")

    # -----------------------------------------------------------------

    def make_young_stars_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making maps of young stars ...")

    # -----------------------------------------------------------------

    def make_ionizing_stars_maps(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Making maps of ionizing stars ...")

# -----------------------------------------------------------------
