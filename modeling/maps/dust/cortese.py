#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.dust.cortese Contains the CorteseDustMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
from astropy import constants

# Import the relevant PTS classes and modules
from ....magic.core.frame import Frame
from ....core.tools import introspection, tables
from ....core.tools import filesystem as fs
from ....core.tools.logging import log
from ..component import MapsComponent
from ....magic.tools.colours import make_colour_map
from ....core.basics.unit import parse_unit as u
from ....core.tools.sequences import combine_unique
from .tirtofuv import TIRtoFUVMapMaker

# -----------------------------------------------------------------

# The path to the table containing the parameters from Cortese et. al 2008
cortese_table_path = fs.join(introspection.pts_dat_dir("modeling"), "cortese.dat")

# The path to the table containing the Galametz calibration parameters
galametz_table_path = fs.join(introspection.pts_dat_dir("modeling"), "galametz.dat")

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

class CorteseDustMapMaker(MapsComponent):

    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(CorteseDustMapMaker, self).__init__(config)

        # -- Attributes --

        # Frames and error maps
        self.frames = dict()
        self.errors = dict()

        # The TIR to FUV map
        self.log_tir_to_fuv = None

        # Maps/dust/cortese path
        self.maps_dust_cortese_path = None

        # The table describing the calibration parameters from Cortese et. al 2008
        # Title of table: Relations to convert the TIR/FUV ratio in A(FUV) for different values of tau and
        # FUV − NIR/optical colours.
        self.cortese = None

        # The SSFR maps (the FUV/optical-NIR colour maps)
        self.ssfr_maps = dict()

        # The attenuation maps (for different FUV/optical-NIR colours)
        self.attenuation_maps = dict()

        # The dust map
        self.map = None

    # -----------------------------------------------------------------

    @classmethod
    def requirements(cls, config=None):

        """
        This function ...
        :param config:
        :return:
        """

        config = cls.get_config(config)

        if config.ssfr_colour == "FUV-H": colour_bands = ["FUV", "2MASS H"]
        elif config.ssfr_colour == "FUV-i": colour_bands = ["FUV", "SDSS i"]
        elif config.ssfr_colour == "FUV-r": colour_bands = ["FUV", "SDSS r"]
        elif config.ssfr_colour == "FUV-g": colour_bands = ["FUV", "SDSS g"]
        elif config.ssfr_colour == "FUV-B": colour_bands = ["FUV", "B"]
        else: raise ValueError("Invalid SSFR colour option")

        # Combine
        return combine_unique(TIRtoFUVMapMaker.requirements(), colour_bands)

    # -----------------------------------------------------------------

    def run(self, log_tir_to_fuv):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup(log_tir_to_fuv)

        # 2. Load the image frames and errors
        self.load_frames()

        # 3. Make the SSFR maps
        self.make_ssfr_maps()

        # 4. Make the dust map
        self.make_map()

        # 5. Make everything positive
        self.make_positive()

        # 6. Normalize the dust map
        self.normalize()

        # 7. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, log_tir_to_fuv):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(CorteseDustMapMaker, self).setup()

        # Load the Cortese et al. 2008 table
        self.cortese = tables.from_file(cortese_table_path, format="ascii.commented_header")

        # Create a maps/dust/cortese directory
        self.maps_dust_cortese_path = fs.create_directory_in(self.maps_dust_path, "cortese")

        # Set the log TIR to FUV map
        self.log_tir_to_fuv = log_tir_to_fuv

    # -----------------------------------------------------------------

    def load_frames(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the necessary data ...")

        data_names = ["GALEX FUV", "2MASS H", "SDSS i", "SDSS r", "SDSS g"]
        #if self.config.ssfr_colour == "FUV-H": data_names.append("2MASS H")
        #elif self.config.ssfr_colour == "FUV-i": data_names.append("SDSS i")
        #elif self.config.ssfr_colour == "FUV-r": data_names.append("SDSS r")
        #elif self.config.ssfr_colour == "FUV-g": data_names.append("SDSS g")
        #elif self.config.ssfr_colour == "FUV-B": data_names.append("")
        #else: raise ValueError("Invalid SSFR colour option")

        # Load all the frames and error maps
        for name in data_names:

            frame = self.dataset.get_frame(name)
            errors = self.dataset.get_errormap(name)

            #from pts.magic.tools import plotting
            #plotting.plot_box(frame)
            #plotting.plot_box(errors)

            # Add the frame and error map to the appropriate dictionary
            self.frames[name] = frame # in original MJy/sr units
            self.errors[name] = errors # in original MJy/sr units

    # -----------------------------------------------------------------

    def make_ssfr_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the sSFR map ...")

        # Loop over the different colour options
        for ssfr_colour in ssfr_colours:

            # Calculate the colour map
            first_band = colour_combinations[ssfr_colour][0]
            second_band = colour_combinations[ssfr_colour][1]
            colour = make_colour_map(self.frames[first_band], self.frames[second_band])

            # Replace NaNs by zeros
            colour.replace_nans(0.0)

            # Mask pixels outside of the low signal-to-noise contour
            #colour[self.mask] = 0.0

            # Set negative pixels to zero
            colour[colour < 0.0] = 0.0

            # Mask low sigal-to-noise pixels in the fuv map, if requested
            #if self.config.ssfr.mask_low_fuv_snr: fuv_h[self.fuv < self.config.ssfr.fuv_snr_level*self.fuv_errors] = 0.0

            # Add the colour map to the dictionary
            self.ssfr_maps[ssfr_colour] = colour

    # -----------------------------------------------------------------

    def make_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the dust map ...")

        # Dust = FUV attenuation = function of (ratio of TIR and FUV luminosity)

        # Loop over the different colour options
        for ssfr_colour in ssfr_colours:

            # Create the FUV attenuation map according to the calibration in Cortese et. al 2008
            fuv_attenuation = self.make_fuv_attenuation_map(ssfr_colour)

            # Set attenuation to zero where the original FUV map is smaller than zero
            fuv_attenuation[self.frames["GALEX FUV"] < 0.0] = 0.0

            # Add the attenuation map to the dictionary
            self.attenuation_maps[ssfr_colour] = fuv_attenuation

        # Choose a specific result as the actual dust map
        self.map = self.attenuation_maps[self.config.ssfr_colour]

    # -----------------------------------------------------------------

    def make_fuv_attenuation_map(self, ssfr_colour):

        """
        This function ...
        :param ssfr_colour
        :return:
        """

        # Inform the user
        log.info("Creating the A(FUV) map according to the relation to the TIR/FUV ratio as described in Cortese et. al 2008 ...")

        ssfr_map = self.ssfr_maps[ssfr_colour]

        # Calculate powers of log(tir_to_fuv)
        tir_to_fuv2 = np.power(self.log_tir_to_fuv, 2.0)
        tir_to_fuv3 = np.power(self.log_tir_to_fuv, 3.0)
        tir_to_fuv4 = np.power(self.log_tir_to_fuv, 4.0)

        # Create an empty image
        a_fuv_cortese = Frame.zeros_like(self.log_tir_to_fuv)

        limits = []
        a1_list = []
        a2_list = []
        a3_list = []
        a4_list = []
        a5_list = []

        # Loop over all entries in the Cortese et. al
        for i in range(len(self.cortese)):

            upper = self.cortese[ssfr_colour][i]
            if i == len(self.cortese) - 1:
                lower = None
            else:
                lower = self.cortese[ssfr_colour][i + 1]

            limits.append((lower, upper))

            a1 = self.cortese["a1"][i]
            a2 = self.cortese["a2"][i]
            a3 = self.cortese["a3"][i]
            a4 = self.cortese["a4"][i]
            a5 = self.cortese["a5"][i]

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

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # 1. Write the SSFR maps
        self.write_ssfr_maps()

        # 2. Write the attenuation maps
        self.write_attenuation_maps()

    # -----------------------------------------------------------------

    def write_ssfr_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the SSFR maps ...")

        # Loop
        for name in self.ssfr_maps:

            # Determine path
            path = fs.join(self.maps_dust_cortese_path, "ssfr_" + name + ".fits")

            # Write
            self.ssfr_maps[name].saveto(path)

    # -----------------------------------------------------------------

    def write_attenuation_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the attenuation maps ...")

        # Loop
        for name in self.attenuation_maps:

            # Determine path
            path = fs.join(self.maps_dust_cortese_path, "fuv_attenuation_" + name + ".fits")

            # Write
            self.attenuation_maps[name].saveto(path)

# -----------------------------------------------------------------
