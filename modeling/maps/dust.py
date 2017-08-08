#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.dust Contains the DustMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import MapsComponent
from ...core.basics.log import log
from ...magic.maps.dust.blackbody import BlackBodyDustMapsMaker
from ...magic.maps.dust.emission import EmissionDustMapsMaker
from ...magic.maps.dust.attenuation import AttenuationDustMapsMaker
from ...magic.maps.dust.hot import HotDustMapsMaker

# -----------------------------------------------------------------

methods = ["black-body", "emission", "attenuation", "hot"]

# -----------------------------------------------------------------

class DustMapMaker(MapsComponent):

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
        super(DustMapMaker, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    @property
    def maps_sub_path(self):

        """
        This function ...
        :return: 
        """

        return self.maps_dust_path

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 3.. Make a dust map based on black body pixel fitting
        if self.config.make_black_body: self.make_black_body()

        # 4. Make a dust map simply based on FIR / submm emission in a certain band
        if self.config.make_emission: self.make_emission()

        # 6. Make a dust map based on UV attenuation
        if self.config.make_attenuation: self.make_attenuation()

        # Make a map of the hot dust
        if self.config.make_hot: self.make_hot()

        # 7. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(DustMapMaker, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def make_black_body(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making a dust map based on black-body fitting to the FIR/submm SED ...")

        # Set method name
        method_name = "black-body"

        # Create the black body dust map maker
        maker = BlackBodyDustMapsMaker(self.config.black_body)

        # Run the maker
        maker.run()

        # Add the dust map to the dictionary
        self.maps[method_name] = maker.maps
        #self.error_maps["black-body"] = maker.error_maps

        # Set origins
        self.origins[method_name] = maker.origins

        # Set methods
        self.methods[method_name] = maker.methods

    # -----------------------------------------------------------------

    def make_emission(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making a dust map based on the emission in a certain FIR/submm band ...")

        # Set method name
        method_name = "emission"

        # Create the emission dust map maker
        maker = EmissionDustMapsMaker()

        # Run the maker
        maker.run()

        # Add the dust map to the dictionary
        self.maps[method_name] = maker.maps

        # Set origins
        self.origins[method_name] = maker.origins

        # Set methods
        self.methods[method_name] = maker.methods

    # -----------------------------------------------------------------

    def make_attenuation(self):

        """
        This function ...
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
        maker.run(attenuation=attenuation_maps, attenuation_origins=attenuation_origins, attenuation_methods=attenuation_methods, method_name=method_name, maps=current)

        # Add the dust maps to the dictionary
        self.maps[method_name] = maker.maps

        # Set origins
        self.origins[method_name] = maker.origins

        # Set methods
        self.methods[method_name] = maker.methods

    # -----------------------------------------------------------------

    def make_hot(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Making a map of the hot dust ...")

        # Set the method name
        method_name = "hot"

        # Get MIPS 24 micron frame and error map
        #mips24 = self.dataset.get_frame("MIPS 24mu")  # in original MJy/sr units
        #mips24_errors = self.dataset.get_errormap("MIPS 24mu")  # in original MJy/sr units
        mips24 = self.get_frame("MIPS 24mu")

        # Get the map of old stars
        #old = self.get_old_stellar_disk_map(self.i1_filter)

        # Get maps of old stars
        old = self.get_old_stellar_disk_maps()
        old_origins = self.get_old_stellar_disk_origins()
        old_methods = self.get_old_stellar_disk_methods()

        # Create the hot dust map maker
        maker = HotDustMapsMaker()

        # Set the factors
        # from 0.2 to 0.7
        factors = self.config.hot_factor_range.linear(self.config.factor_nvalues, as_list=True)

        #print(factors)

        # Get already created maps
        if self.config.remake: current = dict()
        else: current = self.get_current_maps_method(method_name)

        # Run the maker
        #maker.run(mips24=mips24, mips24_errors=mips24_errors, old=old, factors=factors)
        maker.run(mips24=mips24, old=old, old_origins=old_origins, old_methods=old_methods, method_name=method_name, factors=factors, maps=current)

        # Add the dust maps
        self.maps[method_name] = maker.maps

        # Set origins
        self.origins[method_name] = maker.origins

        # Set methods
        self.methods[method_name] = maker.methods

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        #print(self.maps)

        # Write the maps
        self.write_maps()

        # Write origins
        self.write_origins()

        # Write the error maps
        #self.write_error_maps()

        # Write the methods
        self.write_methods()

# -----------------------------------------------------------------
