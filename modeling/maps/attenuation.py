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
from ...magic.maps.attenuation.cortese import CorteseAttenuationMapMaker
from ...magic.maps.attenuation.buat import BuatAttenuationMapMaker

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

        # The maps
        self.maps = dict()

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
        maker = CorteseAttenuationMapMaker()

        # Run the map maker
        maker.run()

        # Set the maps
        self.maps["cortese"] = maker.maps

    # -----------------------------------------------------------------

    def make_buat_attenuation_maps(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Making Buat attenuation maps ...")

        # Create the map maker
        maker = BuatAttenuationMapMaker()

        # Run the map maker
        maker.run()

        # Set the maps
        self.maps["buat"] = maker.maps

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

# -----------------------------------------------------------------
