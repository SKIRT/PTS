#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import MapsAnalysisComponent
from ....core.basics.log import log
from ....magic.maps.dust.attenuation import AttenuationDustMapsMaker
from ....magic.maps.dust.hot import HotDustMapsMaker

# -----------------------------------------------------------------

class DustMapsAnalyser(MapsAnalysisComponent):

    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(DustMapsAnalyser, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        Thisf ucntion ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Make the maps
        self.make_maps()

        # 3. Write
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(DustMapsAnalyser, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This fucntion ...
        :return:
        """

        # Inform the user
        log.info("Making dust maps ...")

        # 1. Make a dust map based on UV attenuation
        self.make_dust_attenuation()

        # 2. Make a map of the hot dust
        self.make_dust_hot()

    # -----------------------------------------------------------------

    def make_dust_attenuation(self):

        """
        Thisnfunction ...
        :return:
        """

        # Inform the user
        log.info("Making a dust map based on the UV attenuation ...")

        # Set the method name
        method_name = "attenuation"

        # Create the Attenuation dust map maker
        maker = AttenuationDustMapsMaker()

        # Get input
        attenuation_maps = self.get_attenuation_maps(flatten=True)
        attenuation_origins = self.get_attenuation_origins(flatten=True)
        attenuation_methods = self.get_attenuation_methods(flatten=True)

        # Get current
        if self.config.remake: current = dict()
        else: current = self.get_current_maps_method(method_name)

        # Run the maker
        maker.run(attenuation=attenuation_maps, attenuation_origins=attenuation_origins,
                  attenuation_methods=attenuation_methods, method_name=method_name, maps=current)

        # Add the dust maps to the dictionary
        self.maps[method_name] = maker.maps

        # Set origins
        self.origins[method_name] = maker.origins

        # Set methods
        self.methods[method_name] = maker.methods

    # -----------------------------------------------------------------

    def make_dust_hot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making a map of the hot dust ...")

        # Set the method name
        method_name = "hot"

        # Get MIPS 24 micron frame and error map
        # mips24 = self.dataset.get_frame("MIPS 24mu")  # in original MJy/sr units
        # mips24_errors = self.dataset.get_errormap("MIPS 24mu")  # in original MJy/sr units
        mips24 = self.get_frame("MIPS 24mu")

        # Get the map of old stars
        # old = self.get_old_stellar_disk_map(self.i1_filter)

        # Get maps of old stars
        old = self.get_old_stellar_disk_maps()
        old_origins = self.get_old_stellar_disk_origins()
        old_methods = self.get_old_stellar_disk_methods()

        # Create the hot dust map maker
        maker = HotDustMapsMaker()

        # Set the factors
        # from 0.2 to 0.7
        factors = self.config.hot_factor_range.linear(self.config.factor_nvalues, as_list=True)

        # print(factors)

        # Get already created maps
        if self.config.remake: current = dict()
        else: current = self.get_current_maps_method(method_name)

        # Run the maker
        # maker.run(mips24=mips24, mips24_errors=mips24_errors, old=old, factors=factors)
        maker.run(mips24=mips24, old=old, old_origins=old_origins, old_methods=old_methods, method_name=method_name,
                  factors=factors, maps=current)

        # Add the dust maps
        self.maps[method_name] = maker.maps

        # Set origins
        self.origins[method_name] = maker.origins

        # Set methods
        self.methods[method_name] = maker.methods

    # -----------------------------------------------------------------

    @property
    def maps_sub_path(self):

        """
        This function ...
        :return:
        """

        return self.dust_path

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the colour maps
        self.write_maps()

        # Write origins
        self.write_origins()

        # Write the methods
        self.write_methods()

# -----------------------------------------------------------------
