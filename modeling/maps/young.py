#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.youngstars Contains the YoungStellarMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from .component import MapsComponent
from ...core.tools import filesystem as fs
from ...core.plot.distribution import DistributionPlotter
from ...magic.core.image import Image
from ...magic.maps.youngstars.young import YoungStellarMapsMaker

# -----------------------------------------------------------------

methods = None

# -----------------------------------------------------------------

class YoungStellarMapMaker(MapsComponent):

    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(YoungStellarMapMaker, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The input FUV and FUV error maps
        self.fuv = None
        self.fuv_errors = None

        # The map of the old stellar disk
        self.old = None

        # The maps of FUV attenuation
        self.fuv_attenuations = None

        # The origins
        self.old_origin = None
        self.fuv_attenuations_origins = None

    # -----------------------------------------------------------------

    @property
    def maps_sub_path(self):

        """
        This function ...
        :return: 
        """

        return self.maps_young_path

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Load the necessary input maps
        self.load_input()

        # 4. Make the map of young stars
        self.make_maps()

        # 5. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(YoungStellarMapMaker, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def load_input(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the necessary input ...")

        # Load the GALEX FUV image and error map
        self.load_fuv()

        # Load FUV attenuation map
        self.load_fuv_attenuation_maps()

        # Load old stellar map
        self.load_old_stellar_map()

    # -----------------------------------------------------------------

    def load_fuv(self):

        """
        This function ...
        :return:
        """

        # Get FUV frame and error map
        self.fuv = self.dataset.get_frame("GALEX FUV") # in original MJy/sr units
        self.fuv_errors = self.dataset.get_errormap("GALEX FUV") # in original MJy/sr units

    # -----------------------------------------------------------------

    def load_fuv_attenuation_maps(self):

        """
        THis function ...
        :return:
        """

        # Inform the user
        log.info("Loading the maps of the FUV attenuation ...")

        # Get the FUV attenuation maps
        self.fuv_attenuations = self.get_fuv_attenuation_maps(flatten=True)

        # Get the FUV attenuation maps origins
        self.fuv_attenuations_origins = self.get_fuv_attenuation_origins(flatten=True)

    # -----------------------------------------------------------------

    def load_old_stellar_map(self):

        """
        This function ...
        """

        # Inform the user
        log.info("Loading the map of old stars ...")

        # Get the map
        self.old = self.get_old_stellar_disk_map(self.i1_filter)

        # Set the old origin
        self.old_origin = self.i1_filter

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Making the maps ...")

        # Create the map maker
        maker = YoungStellarMapsMaker()

        # Set the factors
        factors = self.config.factor_range.linear(self.config.factor_nvalues, as_list=True)

        # Run the map maker
        maker.run(fuv=self.fuv, fuv_errors=self.fuv_errors, old=self.old, fuv_attenuations=self.fuv_attenuations,
                  factors=factors, old_origin=self.old_origin, fuv_attenuations_origins=self.fuv_attenuations_origins)

        # Set the maps
        self.maps = maker.maps

        # Set the origins
        self.origins = maker.origins

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # 1. Write the maps
        self.write_maps()

        # 2. Write origins
        self.write_origins()

# -----------------------------------------------------------------
