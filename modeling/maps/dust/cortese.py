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

# Import the relevant PTS classes and modules
from ....magic.core.frame import Frame
from ....core.tools import introspection, tables
from ....core.tools import filesystem as fs
from ....core.tools.logging import log
from ..component import MapsComponent

# -----------------------------------------------------------------

# The path to the table containing the parameters from Cortese et. al 2008
cortese_table_path = fs.join(introspection.pts_dat_dir("modeling"), "cortese.dat")

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

        # FUV and TIR map in SI units (W/m2)
        self.fuv_si = None
        self.tir_si = None

        # The TIR to FUV ratio map
        self.tir_to_fuv = None

        # The table describing the calibration parameters from Cortese et. al 2008
        # Title of table: Relations to convert the TIR/FUV ratio in A(FUV) for different values of tau and
        # FUV − NIR/optical colours.
        self.cortese = None

        # The output map
        self.map = None

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # Load the image frames and errors
        self.load_frames()

        # Make the FUV map in W/m2 unit
        self.make_fuv()

        # Make the TIR map in W/m2 unit
        self.make_tir()

        # ...
        self.make_map()

        # Normalize the dust map
        self.normalize()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        #ssfr_colour: "FUV-H", "FUV-i", "FUV-r", "FUV-g" or "FUV-B"`
        self.config.ssfr_colour = "FUV-i"

        # Load the Cortese et. al 2008 table
        self.cortese = tables.from_file(cortese_table_path)

    # -----------------------------------------------------------------

    def load_frames(self):

        """
        This function ...
        :return:
        """

        data_names = ["GALEX FUV", "MIPS 24mu", "Pacs blue", "Pacs red"]

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

            self.frames[name] = frame
            self.errors[name] = errors

    # -----------------------------------------------------------------

    def make_fuv(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the FUV map in W/m2 units ...")

        # Convert the FUV map from MJy/sr to W/m2
        exponent = - 20.0 + np.log10(3.e8) - np.log10(0.153e-6) + (2. * np.log10(2.85 / 206264.806247))

        self.fuv_si = self.frames["GALEX FUV"] * 10.0 ** exponent
        self.fuv_si.unit = "W/m2"

        # Save
        #fuv_converted_path = fs.join(self.maps_intermediate_path, "FUV Wpm2.fits")
        #fuv_converted.save(fuv_converted_path)

    # -----------------------------------------------------------------

    def make_tir(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the TIR map in W/m2 units ...")

        # Inform the user
        #log.info("Creating the TIR map in solar units...")

        # MIPS, PACS BLUE AND PACS RED CONVERTED TO LSUN (ABOVE)
        # Galametz (2013) formula for Lsun units
        tir_map = 2.133 * self.images["24mu"].frames.primary + 0.681 * self.images["70mu"].frames.primary + 1.125 * self.images["160mu"].frames.primary
        # Return the TIR map (in solar units)



        # Convert the TIR frame from solar units to W/m2
        exponent = np.log10(3.846e26) - np.log10(4 * np.pi) - (2.0 * np.log10(self.distance_mpc * 3.08567758e22))

        tir_map *= 10.0 ** exponent
        tir_map.unit = Unit("W/m2")

        # Save
        #tir_path = fs.join(self.maps_intermediate_path, "TIR.fits")
        #tir_map.save(tir_path)

    # -----------------------------------------------------------------

    def make_tir_to_fuv(self):

        """
        This function ...
        :return:
        """

        # CALCULATE FUV AND TIR MAP IN W/M2 UNIT

        ## FUV IN W/M2

        ## TIR IN W/M2

        # CALCULATE TIR TO FUV RATIO

        # The ratio of TIR and FUV
        self.tir_to_fuv = np.log10(self.tir_si / self.fuv_si)

        # Save TIR to FUV ratio map
        #tir_to_fuv_path = fs.join(self.maps_intermediate_path, "TIRtoFUV.fits")
        #tir_to_fuv.save(tir_to_fuv_path)

    # -----------------------------------------------------------------

    def make_map(self):

        """
        This function ...
        :return:
        """

        # Dust = FUV attenuation = function of (ratio of TIR and FUV luminosity)

        # Create the FUV attenuation map according to the calibration in Cortese et. al 2008
        a_fuv_cortese = self.create_afuv_cortese("FUV-i")

        # Set attenuation to zero where the original FUV map is smaller than zero
        a_fuv_cortese[self.images["FUV"].frames.primary <= 0.0] = 0.0

        # Make sure all pixel values are larger than or equal to zero
        a_fuv_cortese[a_fuv_cortese < 0.0] = 0.0

        # Cutoff
        a_fuv_cortese[self.cutoff_masks["160mu"]] = 0.0

        # Set the A(FUV) map as the dust map
        self.map = a_fuv_cortese

    # -----------------------------------------------------------------

    def create_afuv_cortese(self, ssfr_colour):

        """
        This function ...
        :param ssfr_colour: "FUV-H", "FUV-i", "FUV-r", "FUV-g" or "FUV-B"
        :return:
        """

        # Inform the user
        log.info("Creating the A(FUV) map according to the relation to the TIR/FUV ratio as described in Cortese et. al 2008 ...")

        # Get the sSFR map
        if ssfr_colour == "FUV-H":
            ssfr = self.get_fuv_h()
        elif ssfr_colour == "FUV-i":
            ssfr = self.get_fuv_i()
        elif ssfr_colour == "FUV-r":
            ssfr = self.get_fuv_r()
        elif ssfr_colour == "FUV-g":
            ssfr = self.get_fuv_g()
        elif ssfr_colour == "FUV-B":
            ssfr = self.get_fuv_b()
        else:
            raise ValueError("Invalid sSFR colour")

        # Calculate powers of tir_to_fuv
        tir_to_fuv2 = np.power(self.tir_to_fuv, 2.0)
        tir_to_fuv3 = np.power(self.tir_to_fuv, 3.0)
        tir_to_fuv4 = np.power(self.tir_to_fuv, 4.0)

        # Create an empty image
        a_fuv_cortese = Frame.zeros_like(self.tir_to_fuv)

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

        # Create the FUV attenuation map
        for i in range(len(limits)):

            if limits[i][0] is None:
                where = ssfr < limits[i][1]
            elif limits[i][1] is None:
                where = ssfr > limits[i][0]
            else:
                where = (ssfr >= limits[i][0]) * (ssfr < limits[i][1])

            # Set the appropriate pixels
            a_fuv_cortese[where] = a1_list[i] + a2_list[i] * self.tir_to_fuv[where] + a3_list[i] * tir_to_fuv2[where] + \
                                   a4_list[i] * tir_to_fuv3[where] + a5_list[i] * tir_to_fuv4[where]

        # Set attenuation to zero where tir_to_fuv is NaN
        a_fuv_cortese[np.isnan(self.tir_to_fuv)] = 0.0

        # Set attenuation to zero where sSFR is smaller than zero
        a_fuv_cortese[ssfr < 0.0] = 0.0

        # Set attenuation to zero where sSFR is greater than 10.5
        a_fuv_cortese[ssfr >= 10.5] = 0.0

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
        fuv_h = -2.5 * np.log10(self.images["FUV"].frames.primary / self.images["H"].frames.primary)

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
        fuv_i = -2.5 * np.log10(self.images["FUV"].frames.primary / self.images["i"].frames.primary)

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
        fuv_r = -2.5 * np.log10(self.images["FUV"].frames.primary / self.images["r"].frames.primary)

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

        pass

    # -----------------------------------------------------------------

    def get_fuv_b(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def normalize(self):

        """
        This function ...
        :return:
        """

        # Normalize the dust map
        self.map.normalize()
        self.map.unit = None

# -----------------------------------------------------------------
