#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.dust.buat Contains the BuatDustMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....core.tools.logging import log
from ..component import MapsComponent
from ....core.tools import filesystem as fs

# -----------------------------------------------------------------

class BuatDustMapMaker(MapsComponent):

    """
    This class...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(BuatDustMapMaker, self).__init__()

        # -- Attributes --

        # Buat parameters
        self.buat = dict()

        # The TIR to FUV ratio map
        self.log_tir_to_fuv = None

        # The dust map
        self.map = None

        # Maps/dust/buat path
        self.maps_dust_buat_path = None

    # -----------------------------------------------------------------

    def run(self, log_tir_to_fuv):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup(log_tir_to_fuv)

        # 2. Make the dust map
        self.make_map()

        # 3. Normalize the dust map
        self.normalize()

        # 4. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, log_tir_to_fuv):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(BuatDustMapMaker, self).setup()

        # Set the buat parameters
        self.buat["NUV"] = (-0.0495, 0.4718, 0.8998, 0.2269)
        self.buat["FUV"] = (-0.0333, 0.3522, 1.1960, 0.4967)

        # Set the TIR to FUV map
        self.log_tir_to_fuv = log_tir_to_fuv

        # Set the path to the maps/dust/buat directory
        self.maps_dust_buat_path = fs.create_directory_in(self.maps_dust_path, "buat")

    # -----------------------------------------------------------------

    def make_map(self):

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

        pass

# -----------------------------------------------------------------
