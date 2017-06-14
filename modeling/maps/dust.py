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
from ...core.tools.logging import log
from ...magic.maps.dust.blackbody import BlackBodyDustMapsMaker
from ...magic.maps.dust.emission import EmissionDustMapsMaker
from ...magic.maps.dust.attenuation import AttenuationDustMapsMaker
from ...magic.maps.dust.hot import HotDustMapsMaker

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

        # Create the black body dust map maker
        maker = BlackBodyDustMapsMaker(self.config.black_body)

        # Run the maker
        maker.run()

        # Add the dust map to the dictionary
        self.maps["black-body"] = maker.maps
        #self.error_maps["black-body"] = maker.error_maps

        # Set origins
        self.origins["black-body"] = maker.origins

    # -----------------------------------------------------------------

    def make_emission(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making a dust map based on the emission in a certain FIR/submm band ...")

        # Create the emission dust map maker
        maker = EmissionDustMapsMaker()

        # Run the maker
        maker.run()

        # Add the dust map to the dictionary
        self.maps["emission"] = maker.maps

        # Set origins
        self.origins["emission"] = maker.origins

    # -----------------------------------------------------------------

    def make_attenuation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making a dust map based on the UV attenuation ...")

        # Create the Attenuation dust map maker
        maker = AttenuationDustMapsMaker()

        # Run the maker
        maker.run()

        # Add the dust maps to the dictionary
        self.maps["attenuation"] = maker.maps

        # Set origins
        self.origins["attenuation"] = maker.origins

    # -----------------------------------------------------------------

    def make_hot(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Making a map of the hot dust ...")

        # Get MIPS 24 micron frame and error map
        mips24 = self.dataset.get_frame("MIPS 24mu")  # in original MJy/sr units
        mips24_errors = self.dataset.get_errormap("MIPS 24mu")  # in original MJy/sr units

        # Get the map of old stars
        old = self.get_old_stellar_disk_map(self.i1_filter)

        # Create the hot dust map maker
        maker = HotDustMapsMaker()

        # Set the factors
        # from 0.2 to 0.7
        factors = self.config.factor_range.linear(self.config.factor_nvalues, as_list=True)

        # Run the maker
        maker.run(mips24=mips24, mips24_errors=mips24_errors, old=old, factors=factors)

        # Add the dust maps
        self.maps["hot"] = maker.maps

        # Set origins
        self.origins["hot"] = maker.origins

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

        # Write origins
        self.write_origins()

        # Write the error maps
        #self.write_error_maps()

# -----------------------------------------------------------------
