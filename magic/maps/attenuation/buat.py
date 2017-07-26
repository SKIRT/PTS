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
from copy import copy

# Import the relevant PTS classes and modules
from ....core.tools.logging import log
from ....core.basics.configurable import Configurable
from ...calibrations.buat import BuatAttenuationCalibration
from .tir_to_uv import make_tir_to_uv
from ...core.frame import Frame
from ....core.filter.filter import parse_filter

# -----------------------------------------------------------------

def make_map(fuv_or_nuv, tir):

    """
    This function ...
    :param fuv_or_nuv:
    :param tir:
    :return:
    """

    # Initialize the map maker
    maker = BuatAttenuationMapsMaker()

    input_dict = dict()
    if fuv_or_nuv.filter == "GALEX FUV": input_dict["fuv"] = fuv_or_nuv
    elif fuv_or_nuv.filter == "GALEX NUV": input_dict["nuv"] = fuv_or_nuv
    else: raise ValueError("First argument must be NUV or FUV map")

    # Set TIR
    input_dict["tirs"] = {"tir": tir}

    # Run the maker
    maker.run(**input_dict)

    # Return the map
    return maker.single_map

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

        # Tirs methods
        self.tirs_methods = None

        # The method name for this class
        self.method_name = None

        # Buat parameters
        self.buat = None

        # The maps
        self.maps = dict()

        # The origins
        self.origins = dict()

        # The methods
        self.methods = dict()

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

        # Origins
        self.tirs_origins = kwargs.pop("tirs_origins", None)

        # Methods
        self.tirs_methods = kwargs.pop("tirs_methods", None)

        # Get the method name for this class
        self.method_name = kwargs.pop("method_name", None)
        if self.has_methods and self.method_name is None: raise ValueError("When methods are specified, method for this class has to be given")

        # Get already created maps
        self.maps = kwargs.pop("maps", dict())

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

    @property
    def has_methods(self):

        """
        This function ...
        :return:
        """

        return self.tirs_methods is not None

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

            # Determine name
            key = "FUV__" + name

            # Set origin
            if self.has_origins:

                # Add FUV as an origin
                origins = copy(self.tirs_origins[name])
                origins.append(parse_filter("FUV"))
                self.origins[key] = origins

            # Set method
            if self.has_methods:

                # Add the method name
                methods = copy(self.tirs_methods[name])
                methods.append(self.method_name)
                self.methods[key] = methods

            # Check whether a map is already present
            if key in self.maps:
                log.warning("The " + name + " attenuation map is already created: not creating it again")
                continue

            # Make the TIR to FUV map
            tir_to_fuv = make_tir_to_uv(self.tirs[name], self.fuv)
            log_tir_to_fuv = Frame(np.log10(tir_to_fuv), wcs=tir_to_fuv.wcs)

            # Calculate FUV attenuation map
            attenuation = parameters[0] * log_tir_to_fuv**3 + parameters[1] * log_tir_to_fuv**2 + parameters[2] * log_tir_to_fuv + parameters[3]

            # Set properties
            attenuation.unit = None # no unit for attenuation
            attenuation.filter = None  # no filter for attenuation
            attenuation.wcs = tir_to_fuv.wcs
            attenuation.distance = tir_to_fuv.distance
            attenuation.pixelscale = tir_to_fuv.pixelscale
            attenuation.psf_filter = tir_to_fuv.psf_filter
            attenuation.fwhm = tir_to_fuv.fwhm

            # Make positive: replace NaNs and negative pixels by zeros
            # Set negatives and NaNs to zero
            attenuation.replace_nans(0.0)
            attenuation.replace_negatives(0.0)

            # Set
            self.maps[key] = attenuation

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

        # Loop over the different TIR maps
        for name in self.tirs:

            # Determine name
            key = "NUV__" + name

            # Set origin
            if self.has_origins:

                # Add NUV as origin
                origins = self.tirs_origins[name]
                origins.append(parse_filter("NUV"))
                self.origins[key] = origins

            # Check whether a map is already present
            if key in self.maps:
                log.warning("The " + name + " attenuation map is already created: not creating it again")
                continue

            # Calculate map
            tir_to_nuv = make_tir_to_uv(self.tirs[name], self.nuv)
            log_tir_to_nuv = Frame(np.log10(tir_to_nuv), wcs=tir_to_nuv.wcs)

            # Calculate attenuation map
            attenuation = parameters[0] * log_tir_to_nuv**3 + parameters[1] * log_tir_to_nuv**2 + parameters[2] * log_tir_to_nuv + parameters[3]

            # Set properties
            attenuation.unit = None # no unit for attenuation
            attenuation.filter = None # no filter for attenuation
            attenuation.wcs = tir_to_nuv.wcs
            attenuation.distance = tir_to_nuv.distance
            attenuation.pixelscale = tir_to_nuv.pixelscale
            attenuation.psf_filter = tir_to_nuv.psf_filter
            attenuation.fwhm = tir_to_nuv.fwhm

            # Make positive: replace NaNs and negative pixels by zeros
            # Set negatives and NaNs to zero
            attenuation.replace_nans(0.0)
            attenuation.replace_negatives(0.0)

            # Set
            self.maps[key] = attenuation

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
