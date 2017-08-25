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
from ...core.basics.log import log
from .component import MapsComponent
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

        # Methods
        self.old_method = None
        self.fuv_attenuations_methods = None

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
        #self.fuv_attenuations, self.fuv_attenuations_origins = self.get_fuv_attenuation_maps_and_origins(flatten=True)
        self.fuv_attenuations, self.fuv_attenuations_origins, self.fuv_attenuations_methods = self.get_fuv_attenuation_maps_origins_and_methods(flatten=True)

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

        # Set the old method
        self.old_method = "disk" #self.get_old_stellar_disk_methods()

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Making the maps ...")

        # Get the current maps
        if self.config.remake: current = dict()
        else: current = self.current_maps

        # Create the map maker
        maker = YoungStellarMapsMaker()

        # Set the factors
        factors = self.config.factor_range.linear(self.config.factor_nvalues, as_list=True)

        # Run the map maker
        maker.run(fuv=self.fuv, fuv_errors=self.fuv_errors, old=self.old, fuv_attenuations=self.fuv_attenuations,
                  factors=factors, old_origin=self.old_origin, fuv_attenuations_origins=self.fuv_attenuations_origins,
                  old_method=self.old_method, fuv_attenuations_methods=self.fuv_attenuations_methods, maps=current)

        # Set the maps
        self.maps = maker.maps

        # Set the origins
        self.origins = maker.origins

        # Set the methods
        self.methods = maker.methods

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

        # Write the methods
        self.write_methods()

# -----------------------------------------------------------------
