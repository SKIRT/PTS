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

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from ....core.tools.logging import log
from ....core.basics.configurable import Configurable
from ...calibrations.buat import BuatAttenuationCalibration
from .tir_to_uv import make_tir_to_uv
from ...core.frame import Frame
from ....core.filter.filter import parse_filter

# -----------------------------------------------------------------

class BuatAttenuationMapsMaker(Configurable):

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
        super(BuatAttenuationMapsMaker, self).__init__(*args, **kwargs)

        # -- Attributes --

        # Input
        self.fuv = None
        self.nuv = None
        self.tirs = None

        # Tirs origins
        self.tirs_origins = None

        # Buat parameters
        self.buat = None

        # The maps
        self.maps = dict()

        # The origins
        self.origins = dict()

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

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(BuatAttenuationMapsMaker, self).setup(**kwargs)

        # Get input
        self.fuv = kwargs.pop("fuv", None)
        self.nuv = kwargs.pop("nuv", None)
        self.tirs = kwargs.pop("tirs")

        # Get distance
        self.distance = kwargs.pop("distance", None)

        # Origins
        self.tirs_origins = kwargs.pop("tirs_origins", None)

        # Create the Cortese instance
        self.buat = BuatAttenuationCalibration()

    # -----------------------------------------------------------------

    @property
    def has_origins(self):

        """
        This function ...
        :return: 
        """

        return self.tirs_origins is not None

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the attenuation maps ...")

        # Make FUV attenuation maps
        if self.fuv is not None: self.make_fuv_maps()

        # Make NUV attenuation maps
        if self.nuv is not None: self.make_nuv_maps()

    # -----------------------------------------------------------------

    def make_fuv_maps(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Creating the FUV attenuation maps ...")

        # Get parameters
        parameters = self.buat.get_fuv_parameters()

        # Loop over the different TIR maps
        for name in self.tirs:

            # Make the TIR to FUV map
            tir_to_fuv = make_tir_to_uv(self.tirs[name], self.fuv, distance=self.distance)
            log_tir_to_fuv = Frame(np.log10(tir_to_fuv), wcs=tir_to_fuv.wcs)

            # Calculate FUV attenuation map
            attenuation = parameters[0] * log_tir_to_fuv**3 + parameters[1] * log_tir_to_fuv**2 + parameters[2] * log_tir_to_fuv + parameters[3]

            # Determine name
            key = name + "_FUV"

            # Make positive: replace NaNs and negative pixels by zeros
            # Set negatives and NaNs to zero
            attenuation.replace_nans(0.0)
            attenuation.replace_negatives(0.0)

            # Set
            self.maps[key] = attenuation

            # Set origin
            if self.has_origins:

                origins = self.tirs_origins[name]
                origins.add(parse_filter("FUV"))
                self.origins[key] = origins

    # -----------------------------------------------------------------

    def make_nuv_maps(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Creating the NUV attenuation maps ...")

        # Get parameters
        parameters = self.buat.get_nuv_parameters()

        # Loop over thed ifferent TIR maps
        for name in self.tirs:

            # Calculate map
            tir_to_nuv = make_tir_to_uv(self.tirs[name], self.nuv, distance=self.distance)
            log_tir_to_nuv = Frame(np.log10(tir_to_nuv), wcs=tir_to_nuv.wcs)

            # Calculate attenuation map
            attenuation = parameters[0] * log_tir_to_nuv**3 + parameters[1] * log_tir_to_nuv**2 + parameters[2] * log_tir_to_nuv + parameters[3]

            # Determine name
            key = name + "_NUV"

            # Make positive: replace NaNs and negative pixels by zeros
            # Set negatives and NaNs to zero
            attenuation.replace_nans(0.0)
            attenuation.replace_negatives(0.0)

            # Set
            self.maps[key] = attenuation

            # Set origin
            if self.has_origins:

                origins = self.tirs_origins[name]
                origins.add(parse_filter("NUV"))
                self.origins[key] = origins

    # -----------------------------------------------------------------

    @property
    def single_map(self):

        """
        This function ...
        :return: 
        """

        if len(self.maps) != 1: raise ValueError("Not a single map")
        return self.maps[self.maps.keys()[0]]

# -----------------------------------------------------------------
