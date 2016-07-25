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
from astropy.units import Unit
from astropy import constants

# Import the relevant PTS classes and modules
from ....magic.core.frame import Frame
from ....core.tools import introspection, tables
from ....core.tools import filesystem as fs
from ....core.tools.logging import log
from ..component import MapsComponent

# -----------------------------------------------------------------

# The path to the table containing the parameters from Cortese et. al 2008
cortese_table_path = fs.join(introspection.pts_dat_dir("modeling"), "cortese.dat")

# The path to the table containing the Galametz calibration parameters
galametz_table_path = fs.join(introspection.pts_dat_dir("modeling"), "galametz.dat")

# -----------------------------------------------------------------

speed_of_light = constants.c
solar_luminosity = 3.846e26 * Unit("W")

# -----------------------------------------------------------------

class CorteseDustMapMaker(MapsComponent):

    """
    This class...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(CorteseDustMapMaker, self).__init__()

        # -- Attributes --

        # Frames and error maps
        self.frames = dict()
        self.errors = dict()

        # The TIR to FUV map
        self.log_tir_to_fuv = None

        # The sSFR map
        self.ssfr = None

        # Maps/dust/cortese path
        self.maps_dust_cortese_path = None

        # The table describing the calibration parameters from Cortese et. al 2008
        # Title of table: Relations to convert the TIR/FUV ratio in A(FUV) for different values of tau and
        # FUV − NIR/optical colours.
        self.cortese = None

        # The dust map
        self.map = None

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

        # 3. Make the sSFR map
        self.make_ssfr()

        # 4. Make the dust map
        self.make_map()

        # 5. Normalize the dust map
        self.normalize()

        # 6. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, log_tir_to_fuv):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(CorteseDustMapMaker, self).setup()

        #ssfr_colour: "FUV-H", "FUV-i", "FUV-r", "FUV-g" or "FUV-B"`
        #self.config.ssfr_colour = "FUV-i"
        self.config.ssfr_colour = "FUV-H"

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

        data_names = ["GALEX FUV"]
        if self.config.ssfr_colour == "FUV-H": data_names.append("2MASS H")
        elif self.config.ssfr_colour == "FUV-i": data_names.append("SDSS i")
        elif self.config.ssfr_colour == "FUV-r": data_names.append("SDSS r")
        elif self.config.ssfr_colour == "FUV-g": data_names.append("SDSS g")
        #elif self.config.ssfr_colour == "FUV-B": data_names.append("")
        else: raise ValueError("Invalid SSFR colour option")

        # Load all the frames and error maps
        for name in data_names:

            frame = self.dataset.get_frame(name)
            errors = self.dataset.get_errors(name)

            # Add the frame and error map to the appropriate dictionary
            self.frames[name] = frame # in original MJy/sr units
            self.errors[name] = errors # in original MJy/sr units

    # -----------------------------------------------------------------

    def make_ssfr(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the sSFR map ...")

        # Get the sSFR map
        if self.config.ssfr_colour == "FUV-H": self.ssfr = self.get_fuv_h()
        elif self.config.ssfr_colour == "FUV-i": self.ssfr = self.get_fuv_i()
        elif self.config.ssfr_colour == "FUV-r": self.ssfr = self.get_fuv_r()
        elif self.config.ssfr_colour == "FUV-g": self.ssfr = self.get_fuv_g()
        #elif self.config.ssfr_colour == "FUV-B": self.ssfr = self.get_fuv_b()
        else: raise ValueError("Invalid sSFR colour")

    # -----------------------------------------------------------------

    def make_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the dust map ...")

        # Dust = FUV attenuation = function of (ratio of TIR and FUV luminosity)

        # Create the FUV attenuation map according to the calibration in Cortese et. al 2008
        a_fuv_cortese = self.create_afuv_cortese()

        # Set attenuation to zero where the original FUV map is smaller than zero
        a_fuv_cortese[self.frames["GALEX FUV"] <= 0.0] = 0.0

        # Make sure all pixel values are larger than or equal to zero
        #a_fuv_cortese[a_fuv_cortese < 0.0] = 0.0

        # Cutoff
        #a_fuv_cortese[self.cutoff_masks["160mu"]] = 0.0

        # Set the A(FUV) map as the dust map
        self.map = a_fuv_cortese

    # -----------------------------------------------------------------

    def create_afuv_cortese(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the A(FUV) map according to the relation to the TIR/FUV ratio as described in Cortese et. al 2008 ...")

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

            upper = self.cortese[self.config.ssfr_colour][i]
            if i == len(self.cortese) - 1:
                lower = None
            else:
                lower = self.cortese[self.config.ssfr_colour][i + 1]

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
                where = self.ssfr < upper_limit
            elif upper_limit is None:
                where = self.ssfr > lower_limit
            else: where = (self.ssfr >= lower_limit) * (self.ssfr < upper_limit)

            # Set the appropriate pixels
            a_fuv_cortese[where] = a1_list[i] + a2_list[i] * self.log_tir_to_fuv[where] + a3_list[i] * tir_to_fuv2[where] + \
                                   a4_list[i] * tir_to_fuv3[where] + a5_list[i] * tir_to_fuv4[where]

        # The absolute upper limit (so 10.5 for FUV-H, 7.5 for FUV-i, 7.3 for FUV-r, 6.7 for FUV-g, and 6.3 for FUV-B
        absolute_upper_limit = limits[0][1]

        # Set attenuation to zero where tir_to_fuv is NaN
        a_fuv_cortese[np.isnan(self.log_tir_to_fuv)] = 0.0

        # Set attenuation to zero where sSFR is smaller than zero
        a_fuv_cortese[self.ssfr < 0.0] = 0.0

        # Set attenuation to zero where sSFR is greater than the absolute upper limit for the FUV-IR/optical colour
        a_fuv_cortese[self.ssfr >= absolute_upper_limit] = 0.0

        # Return the A(FUV) map
        return a_fuv_cortese

    # -----------------------------------------------------------------

    def get_fuv_h(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the FUV-H colour map ...")

        # Calculate the colour map
        fuv_h = Frame(-2.5 * np.log10(self.frames["GALEX FUV"] / self.frames["2MASS H"]))

        # Replace NaNs by zeros
        fuv_h.replace_nans(0.0)

        # Mask pixels outside of the low signal-to-noise contour
        # fuv_h[self.mask] = 0.0

        # Set negative pixels to zero
        fuv_h[fuv_h < 0.0] = 0.0

        # Mask low sigal-to-noise pixels in the fuv map, if requested
        # if self.config.ssfr.mask_low_fuv_snr: fuv_h[self.fuv < self.config.ssfr.fuv_snr_level*self.fuv_errors] = 0.0

        # Return the FUV-H colour map
        return fuv_h

    # -----------------------------------------------------------------

    def get_fuv_i(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the FUV-i colour map ...")

        # Calculate the colour map
        fuv_i = Frame(-2.5 * np.log10(self.frames["GALEX FUV"] / self.frames["SDSS i"]))

        # Replace NaNs by zeros
        fuv_i.replace_nans(0.0)

        # Mask pixels outside of the low signal-to-noise contour
        # fuv_i[self.mask] = 0.0

        # Set negative pixels to zero
        fuv_i[fuv_i < 0.0] = 0.0

        # Return the FUV-i colour map
        return fuv_i

    # -----------------------------------------------------------------

    def get_fuv_r(self):

        """
        This function ...
        :return:
        """

        # Calculate the colour map
        fuv_r = Frame(-2.5 * np.log10(self.frames["GALEX FUV"] / self.frames["SDSS r"]))

        # Replace NaNs by zeros
        fuv_r.replace_nans(0.0)

        # Mask pixels outside of the low signal-to-noise contour
        # fuv_r[self.mask] = 0.0

        # Set negative pixels to zero
        fuv_r[fuv_r < 0.0] = 0.0

        # Return the FUV-r colour map
        return fuv_r

    # -----------------------------------------------------------------

    def get_fuv_g(self):

        """
        This function ...
        :return:
        """

        # Calculate the colour map
        fuv_g = Frame(-2.5 * np.log10(self.frames["GALEX FUV"] / self.frames["SDSS g"]))

        # Replace NaNs by zeros
        fuv_g.replace_nans(0.0)

        # Mask pixels outside of the low signal-to-noise contour
        # fuv_r[self.mask] = 0.0

        # Set negative pixels to zero
        fuv_g[fuv_g < 0.0] = 0.0

        # Return the FUV-g colour map
        return fuv_g

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

        # Write the sSFR map
        self.write_ssfr()

    # -----------------------------------------------------------------

    def write_ssfr(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the sSFR map ...")

        # Write the sSFR map
        path = fs.join(self.maps_dust_cortese_path, "ssfr.fits")

        # Write
        self.ssfr.save(path)

# -----------------------------------------------------------------
