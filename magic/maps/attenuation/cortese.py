#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.maps.dust.attenuation Contains the AttenuationDustMapsMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from ....magic.core.frame import Frame
from ....core.tools.logging import log
from ....core.basics.configurable import Configurable
from ...calibrations.cortese import CorteseAttenuationCalibration
from .tir_to_uv import make_tir_to_uv
from ....core.filter.filter import parse_filter
from ....core.tools import sequences

# -----------------------------------------------------------------

def make_map(fuv, tir, ssfr, ssfr_colour):

    """
    This function ...
    :param fuv:
    :param tir:
    :param ssfr:
    :param ssfr_colour:
    :return: 
    """

    # Create the attenuation map maker
    maker = CorteseAttenuationMapsMaker()

    # Set input
    tirs = {"standard": tir}
    ssfrs = {ssfr_colour: ssfr}

    # Run
    maker.run(fuv=fuv, tirs=tirs, ssfrs=ssfrs)

    # Get the map
    return maker.single_map

# -----------------------------------------------------------------

class CorteseAttenuationMapsMaker(Configurable):

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
        super(CorteseAttenuationMapsMaker, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The FUV map
        self.fuv = None

        # The TIR maps
        self.tirs = None

        # The ssfr maps
        self.ssfrs = None

        # Origins
        self.tirs_origins = None
        self.ssfrs_origins = None

        # The table describing the calibration parameters from Cortese et. al 2008
        # Title of table: Relations to convert the TIR/FUV ratio in A(FUV) for different values of tau and
        # FUV − NIR/optical colours.
        self.cortese = None

        # The SSFR maps (the FUV/optical-NIR colour maps)
        self.ssfrs = dict()

        # The attenuation maps (for different FUV/optical-NIR colours)
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
        super(CorteseAttenuationMapsMaker, self).setup(**kwargs)

        # Get input
        self.fuv = kwargs.pop("fuv")
        self.tirs = kwargs.pop("tirs")
        self.ssfrs = kwargs.pop("ssfrs")

        # Get origins
        self.tirs_origins = kwargs.pop("tirs_origins", None)
        self.ssfrs_origins = kwargs.pop("ssfrs_origins", None)

        # Create the Cortese instance
        self.cortese = CorteseAttenuationCalibration()

    # -----------------------------------------------------------------

    @property
    def has_origins(self):

        """
        This function ...
        :return: 
        """

        return self.tirs_origins is not None and self.ssfrs_origins is not None

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the attenuation maps ...")

        # Dust = FUV attenuation = function of (ratio of TIR and FUV luminosity)

        # Loop over the different TIR maps
        for name in self.tirs:

            # Make the TIR to FUV map
            tir_to_fuv = make_tir_to_uv(self.tirs[name], self.fuv)
            log_tir_to_fuv = Frame(np.log10(tir_to_fuv), wcs=tir_to_fuv.wcs)

            # Loop over the different colour options
            for ssfr_colour in self.ssfrs:

                # Get the ssfr map
                ssfr = self.ssfrs[ssfr_colour]

                # Create the FUV attenuation map according to the calibration in Cortese et. al 2008
                fuv_attenuation = make_fuv_attenuation_map(self.cortese, ssfr_colour, log_tir_to_fuv, ssfr)

                # Set attenuation to zero where the original FUV map is smaller than zero
                fuv_attenuation[self.fuv < 0.0] = 0.0

                # Determine name
                key = name + "_" + ssfr_colour

                # Make positive: replace NaNs and negative pixels by zeros
                # Set negatives and NaNs to zero
                fuv_attenuation.replace_nans(0.0)
                fuv_attenuation.replace_negatives(0.0)

                # Add the attenuation map to the dictionary
                self.maps[key] = fuv_attenuation

                # Set origins
                if self.has_origins:

                    origins = self.tirs_origins[name]
                    sequences.extend_unique(origins, self.ssfrs_origins[ssfr_colour])
                    sequences.append_unique(origins, parse_filter("FUV"))
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

def make_fuv_attenuation_map(cortese, ssfr_colour, log_tir_to_fuv, ssfr):

    """
    This function ...
    :param cortese:
    :param ssfr_colour:
    :param log_tir_to_fuv:
    :param ssfr:
    :return:
    """

    # Inform the user
    log.info("Creating the A(FUV) map according to the relation to the TIR/FUV ratio as described in Cortese et. al 2008 ...")

    # Calculate powers of log(tir_to_fuv)
    tir_to_fuv2 = np.power(log_tir_to_fuv.data, 2.0)
    tir_to_fuv3 = np.power(log_tir_to_fuv.data, 3.0)
    tir_to_fuv4 = np.power(log_tir_to_fuv.data, 4.0)

    # Create an empty image
    a_fuv_cortese = Frame.zeros_like(log_tir_to_fuv)

    # Create the FUV attenuation map
    for tau, colour_range, parameters in cortese.taus_ranges_and_parameters(ssfr_colour):

        # Debugging
        log.debug("Setting FUV attenuation values for tau = " + str(tau) + " ...")

        # Set mask
        where = (ssfr >= colour_range.min) * (ssfr < colour_range.max)

        # Set the appropriate pixels
        a_fuv_cortese[where] = parameters[0] + parameters[1] * log_tir_to_fuv[where] + parameters[2] * tir_to_fuv2[where] + \
                               parameters[3] * tir_to_fuv3[where] + parameters[4] * tir_to_fuv4[where]

    # Get absolute upper limit
    absolute_upper_limit = cortese.upper_limit(ssfr_colour)

    # Set attenuation to zero where tir_to_fuv is NaN
    a_fuv_cortese[np.isnan(log_tir_to_fuv)] = 0.0

    # Set attenuation to zero where sSFR is smaller than zero
    a_fuv_cortese[ssfr < 0.0] = 0.0

    # Set attenuation to zero where sSFR is greater than the absolute upper limit for the FUV-IR/optical colour
    a_fuv_cortese[ssfr >= absolute_upper_limit] = 0.0

    # Return the A(FUV) map
    return a_fuv_cortese

# -----------------------------------------------------------------
