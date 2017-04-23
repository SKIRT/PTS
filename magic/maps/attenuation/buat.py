#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.dust.buat Contains the BuatDustMapsMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....core.tools.logging import log
from ....core.basics.configurable import Configurable

# -----------------------------------------------------------------

class BuatAttenuationMapsMaker(Configurable):

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
        super(BuatAttenuationMapsMaker, self).__init__(config, interactive)

        # -- Attributes --

        # Buat parameters
        self.buat = dict()

        # The TIR to FUV ratio map
        self.log_tir_to_fuv = None

        # The dust map
        #self.map = None

        # The maps
        self.maps = dict()

        # Maps/dust/buat path
        self.maps_dust_buat_path = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Make the dust map
        self.make_maps()

        # 3. Make everything positive
        #self.make_positive()

        # 4. Normalize the dust map
        #self.normalize()

        # 5. Writing
        #self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(BuatAttenuationMapsMaker, self).setup()

        # Set the buat parameters
        self.buat["NUV"] = (-0.0495, 0.4718, 0.8998, 0.2269)
        self.buat["FUV"] = (-0.0333, 0.3522, 1.1960, 0.4967)

        # Set the TIR to FUV map
        self.log_tir_to_fuv = log_tir_to_fuv

        # Set the path to the maps/dust/buat directory
        #self.maps_dust_buat_path = fs.create_directory_in(self.maps_dust_path, "buat")

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the dust map ...")

        # Calculate FUV attenuation map
        a_fuv = -0.0333 * self.log_tir_to_fuv**3 + 0.3522 * self.log_tir_to_fuv**2 + 1.1960 * self.log_tir_to_fuv + 0.4967

        # Set map
        self.map = a_fuv

    # -----------------------------------------------------------------

    def make_positive(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Replacing NaNs and negative pixels by zeros ...")

        # Set negatives and NaNs to zero
        self.map.replace_nans(0.0)
        self.map.replace_negatives(0.0)

    # -----------------------------------------------------------------

    def normalize(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Normalizing the dust map ...")

        # Normalize the dust map
        self.map.normalize()
        self.map.unit = None

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

# -----------------------------------------------------------------
