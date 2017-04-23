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

# Import astronomical modules
from astropy import constants

# Import the relevant PTS classes and modules
from ....magic.core.frame import Frame
from ....core.tools import filesystem as fs
from ....core.tools.logging import log
from ....magic.tools.colours import make_colour_map
from ....core.units.parsing import parse_unit as u
from ....core.basics.configurable import Configurable

# -----------------------------------------------------------------

speed_of_light = constants.c
solar_luminosity = 3.846e26 * u("W")

# -----------------------------------------------------------------

#ssfr_colours = ["FUV-H", "FUV-i", "FUV-r", "FUV-g", "FUV-B"]
ssfr_colours = ["FUV-H", "FUV-i", "FUV-r", "FUV-g"]

colour_combinations = {"FUV-H": ("GALEX FUV", "2MASS H"),
                       "FUV-i": ("GALEX FUV", "SDSS i"),
                       "FUV-r": ("GALEX FUV", "SDSS r"),
                       "FUV-g": ("GALEX FUV", "SDSS g")}

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

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(CorteseAttenuationMapsMaker, self).__init__(config, interactive)

        # -- Attributes --

        # The FUV map
        self.fuv = None

        # The TIR maps
        self.tirs = None

        # The ssfr maps
        self.ssfrs = None

        # Frames and error maps
        #self.frames = dict()
        #self.errors = dict()

        # The TIR to FUV map
        #self.log_tir_to_fuv = None

        # Maps/dust/cortese path
        #self.maps_dust_cortese_path = None

        # The table describing the calibration parameters from Cortese et. al 2008
        # Title of table: Relations to convert the TIR/FUV ratio in A(FUV) for different values of tau and
        # FUV − NIR/optical colours.
        self.cortese = None

        # The SSFR maps (the FUV/optical-NIR colour maps)
        self.ssfrs = dict()

        # The attenuation maps (for different FUV/optical-NIR colours)
        self.maps = dict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Load the image frames and errors
        #self.load_frames()

        # 3. Make the SSFR maps
        #self.make_ssfr_maps()

        # 4. Make the dust map
        self.make_maps()

        # 5. Make everything positive
        #self.make_positive()

        # 6. Normalize the dust map
        #self.normalize()

        # 7. Writing
        #self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(CorteseAttenuationMapsMaker, self).setup(**kwargs)

        # Load the Cortese et al. 2008 table
        #self.cortese = tables.from_file(cortese_table_path, format="ascii.commented_header")

        # Create a maps/dust/cortese directory
        #self.maps_dust_cortese_path = fs.create_directory_in(self.maps_dust_path, "cortese")

        # Set the log TIR to FUV map
        #self.log_tir_to_fuv = log_tir_to_fuv

        # Get input
        self.fuv = kwargs.pop("fuv")
        self.tirs = kwargs.pop("tirs")
        self.ssfrs = kwargs.pop("ssfrs")

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the dust map ...")

        # Dust = FUV attenuation = function of (ratio of TIR and FUV luminosity)

        # Loop over the different TIR maps
        for name in self.tirs:

            # Make the TIR to FUV map
            tir_to_fuv = make_tir_to_fuv_map(self.tirs[name], self.fuv)
            log_tir_to_fuv = Frame(np.log10(tir_to_fuv))

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

                # Add the attenuation map to the dictionary
                self.maps[key] = fuv_attenuation

        # Choose a specific result as the actual dust map
        #self.map = self.attenuation_maps[self.config.ssfr_colour]

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

    def make_positive(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Replacing NaNs and negative pixels by zeros ...")

        # Set negatives and NaNs to zero
        self.map.replace_nans(0.0)
        self.map.replace_negatives(0.0)

    # -----------------------------------------------------------------

    def normalize(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Normalizing the dust map ...")

        # Normalize the dust map
        self.map.normalize()
        self.map.unit = None

# -----------------------------------------------------------------

def make_tir_to_fuv_map(tir, fuv):

    """
    This function ...
    :param tir: 
    :param fuv: 
    :return: 
    """

    # Conversions necessary? -> YES!

    ## Convert the FUV map from Lsun to W/m2
    #assert self.frames["GALEX FUV"].unit == "Lsun"
    ## Convert the TIR map from Lsun to W / m2
    #conversion_factor = 1.0
    # Conversion from Lsun to W
    #conversion_factor *= solar_luminosity.to("W").value
    # Conversion from W [LUMINOSITY] to W / m2 [FLUX]
    #distance = self.galaxy_properties.distance
    #conversion_factor /= (4. * np.pi * distance ** 2).to("m2").value
    # FUV in W/M2
    #self.fuv_si = self.frames["GALEX FUV"] * conversion_factor
    #self.fuv_si.unit = "W/m2"

    ## Convert the TIR map from Lsun to W / m2
    #conversion_factor = 1.0
    # Conversion from Lsun to W
    #conversion_factor *= solar_luminosity.to("W").value
    # Conversion from W [LUMINOSITY] to W / m2 [FLUX]
    #distance = self.galaxy_properties.distance
    #conversion_factor /= (4. * np.pi * distance ** 2).to("m2").value
    ## CONVERT AND SET NEW UNIT
    #self.tir_si = Frame(tir_map * conversion_factor)
    #self.tir_si.unit = "W/m2"

    # CALCULATE FUV AND TIR MAP IN W/M2 UNIT

    ## FUV IN W/M2

    ## TIR IN W/M2

    # CALCULATE TIR TO FUV RATIO

    # The ratio of TIR and FUV
    tir_to_fuv = tir_si / fuv_si
    #log_tir_to_fuv = Frame(np.log10(self.tir_to_fuv)

    return tir_to_fuv

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

    #ssfr_map = self.ssfr_maps[ssfr_colour]

    # Calculate powers of log(tir_to_fuv)
    tir_to_fuv2 = np.power(log_tir_to_fuv, 2.0)
    tir_to_fuv3 = np.power(log_tir_to_fuv, 3.0)
    tir_to_fuv4 = np.power(log_tir_to_fuv, 4.0)

    # Create an empty image
    a_fuv_cortese = Frame.zeros_like(log_tir_to_fuv)

    limits = []
    a1_list = []
    a2_list = []
    a3_list = []
    a4_list = []
    a5_list = []

    # Loop over all entries in the Cortese et. al
    for i in range(len(cortese)):

        upper = cortese[ssfr_colour][i]
        if i == len(cortese) - 1: lower = None
        else: lower = cortese[ssfr_colour][i + 1]

        limits.append((lower, upper))

        a1 = cortese["a1"][i]
        a2 = cortese["a2"][i]
        a3 = cortese["a3"][i]
        a4 = cortese["a4"][i]
        a5 = cortese["a5"][i]

        a1_list.append(a1)
        a2_list.append(a2)
        a3_list.append(a3)
        a4_list.append(a4)
        a5_list.append(a5)

    # Debugging
    log.debug("a1 values: " + " ".join([str(a) for a in a1_list]))
    log.debug("a2 values: " + " ".join([str(a) for a in a2_list]))
    log.debug("a3 values: " + " ".join([str(a) for a in a3_list]))
    log.debug("a4 values: " + " ".join([str(a) for a in a4_list]))
    log.debug("a5 values: " + " ".join([str(a) for a in a5_list]))

    # Create the FUV attenuation map
    for i in range(len(limits)):

        upper_limit = limits[i][1]
        lower_limit = limits[i][0]

        if lower_limit is None:
            where = ssfr_map < upper_limit
        elif upper_limit is None:
            where = ssfr_map > lower_limit
        else: where = (ssfr_map >= lower_limit) * (ssfr_map < upper_limit)

        # Set the appropriate pixels
        a_fuv_cortese[where] = a1_list[i] + a2_list[i] * self.log_tir_to_fuv[where] + a3_list[i] * tir_to_fuv2[where] + \
                               a4_list[i] * tir_to_fuv3[where] + a5_list[i] * tir_to_fuv4[where]

    # The absolute upper limit (so 10.5 for FUV-H, 7.5 for FUV-i, 7.3 for FUV-r, 6.7 for FUV-g, and 6.3 for FUV-B
    absolute_upper_limit = limits[0][1]

    # Set attenuation to zero where tir_to_fuv is NaN
    a_fuv_cortese[np.isnan(self.log_tir_to_fuv)] = 0.0

    # Set attenuation to zero where sSFR is smaller than zero
    a_fuv_cortese[ssfr_map < 0.0] = 0.0

    # Set attenuation to zero where sSFR is greater than the absolute upper limit for the FUV-IR/optical colour
    a_fuv_cortese[ssfr_map >= absolute_upper_limit] = 0.0

    # Return the A(FUV) map
    return a_fuv_cortese

# -----------------------------------------------------------------
