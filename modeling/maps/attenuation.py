#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.attenuation Contains the AttenuationMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from .component import MapsComponent
from ...magic.maps.attenuation.cortese import CorteseAttenuationMapsMaker
from ...magic.maps.attenuation.buat import BuatAttenuationMapsMaker

# -----------------------------------------------------------------

class AttenuationMapMaker(MapsComponent):

    """
    This class...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(AttenuationMapMaker, self).__init__(config, interactive)

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Setup
        self.setup(**kwargs)

        # 2. Cortese
        self.make_cortese_attenuation_maps()

        # 3. Buat
        self.make_buat_attenuation_maps()

        # 4. Write
        self.write()

    # -----------------------------------------------------------------

    def make_cortese_attenuation_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making Cortese attenuation maps ...")

        # Create the map maker
        maker = CorteseAttenuationMapsMaker()

        # Get the input
        fuv = self.get_frame_for_filter(self.fuv_filter)
        tirs = self.get_tir_maps()
        ssfrs = self.get_ssfr_maps()
        tirs_origins = self.get_tir_origins()
        ssfrs_origins = self.get_ssfr_origins()

        # Run the map maker
        maker.run(fuv=fuv, tirs=tirs, ssfrs=ssfrs, tirs_origins=tirs_origins, ssfrs_origins=ssfrs_origins)

        # Set the maps
        self.maps["cortese"] = maker.maps

        # Set the origins
        self.maps["cortese"] = maker.origins

    # -----------------------------------------------------------------

    def make_buat_attenuation_maps(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Making Buat attenuation maps ...")

        # Create the map maker
        maker = BuatAttenuationMapsMaker()

        # Get the input
        fuv = self.get_frame_for_filter(self.fuv_filter)
        nuv = self.get_frame_for_filter(self.nuv_filter)
        tirs = self.get_tir_maps()
        tirs_origins = self.get_tir_origins()

        # Run the map maker
        maker.run(fuv=fuv, nuv=nuv, tirs=tirs, tirs_origins=tirs_origins)

        # Set the maps
        self.maps["buat"] = maker.maps

        # Set the origins
        self.maps["origins"] = maker.origins

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the maps
        self.write_maps()

        # Write the origins
        self.write_origins()

# -----------------------------------------------------------------
