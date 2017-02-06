#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.dust.tirtofuv Contains the TIRtoFUVMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
from astropy import constants

# Import the relevant PTS classes and modules
from ....core.tools import introspection, tables
from ....core.tools import filesystem as fs
from ....core.tools.logging import log
from ..component import MapsComponent
from ....magic.core.frame import Frame
from ....core.basics.unit import parse_unit as u

# -----------------------------------------------------------------

# The path to the table containing the parameters from Cortese et. al 2008
cortese_table_path = fs.join(introspection.pts_dat_dir("modeling"), "cortese.dat")

# The path to the table containing the Galametz calibration parameters
galametz_table_path = fs.join(introspection.pts_dat_dir("modeling"), "galametz.dat")

# -----------------------------------------------------------------

speed_of_light = constants.c
solar_luminosity = 3.846e26 * u("W")

# -----------------------------------------------------------------

class TIRtoFUVMapMaker(MapsComponent):

    """
    This class...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(TIRtoFUVMapMaker, self).__init__()

        # -- Attributes --

        # Maps/dust/tirfuv path
        self.maps_tirfuv_path = None

        # Frames and error maps
        self.frames = dict()
        self.errors = dict()

        # FUV and TIR map in SI units (W/m2)
        self.fuv_si = None
        self.tir_si = None

        # The TIR to FUV ratio map
        self.tir_to_fuv = None
        self.log_tir_to_fuv = None

        # The table with the calibration factors from galametz (TIR ..)
        self.galametz = None

    # -----------------------------------------------------------------

    @classmethod
    def requirements(cls, config=None):

        """
        This function ...
        :param config:
        :return:
        """

        config = cls.get_config(config)
        return ["GALEX FUV", "MIPS 24mu", "Pacs blue", "Pacs red"]

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the image frames and errors
        self.load_frames()

        # 3. Make the FUV map in W/m2 unit
        self.make_fuv()

        # 4. Make the TIR map in W/m2 unit
        self.make_tir()

        # 5. Make the TIR to FUV ratio map
        self.make_tir_to_fuv()

        # 6. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(TIRtoFUVMapMaker, self).setup()

        # Load the Galametz et al. table
        self.galametz = tables.from_file(galametz_table_path, format="ascii.commented_header")

        # Create a maps/dust/cortese directory
        self.maps_tirfuv_path = fs.create_directory_in(self.maps_dust_path, "tir fuv")

    # -----------------------------------------------------------------

    def load_frames(self):

        """
        This function ...
        :return:
        """

        # Get the galaxy distance
        distance = self.galaxy_properties.distance

        # Load all the frames and error maps
        for name in self.requirements(self.config):

            frame = self.dataset.get_frame(name)
            errors = self.dataset.get_errormap(name)

            ## CONVERT TO LSUN

            # Get pixelscale and wavelength
            pixelscale = frame.average_pixelscale
            wavelength = frame.filter.pivot

            ##

            # Conversion from MJy / sr to Jy / sr
            conversion_factor = 1e6

            # Conversion from Jy / sr to Jy / pix(2)
            conversion_factor *= (pixelscale ** 2).to("sr/pix2").value

            # Conversion from Jy / pix to W / (m2 * Hz) (per pixel)
            conversion_factor *= 1e-26

            # Conversion from W / (m2 * Hz) (per pixel) to W / (m2 * m) (per pixel)
            conversion_factor *= (speed_of_light / wavelength**2).to("Hz/m").value

            # Conversion from W / (m2 * m) (per pixel) [SPECTRAL FLUX] to W / m [SPECTRAL LUMINOSITY]
            conversion_factor *= (4. * np.pi * distance**2).to("m2").value

            # Conversion from W / m [SPECTRAL LUMINOSITY] to W [LUMINOSITY]
            conversion_factor *= wavelength.to("m").value

            # Conversion from W to Lsun
            conversion_factor *= 1. / solar_luminosity.to("W").value

            ## CONVERT

            frame *= conversion_factor
            frame.unit = "Lsun"

            errors *= conversion_factor
            errors.unit = "Lsun"

            # Add the frame and error map to the appropriate dictionary
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

        ## Convert the FUV map from Lsun to W/m2

        assert self.frames["GALEX FUV"].unit == "Lsun"

        ## Convert the TIR map from Lsun to W / m2

        conversion_factor = 1.0

        # Conversion from Lsun to W

        conversion_factor *= solar_luminosity.to("W").value

        # Conversion from W [LUMINOSITY] to W / m2 [FLUX]
        distance = self.galaxy_properties.distance
        conversion_factor /= (4. * np.pi * distance ** 2).to("m2").value

        # FUV in W/M2
        self.fuv_si = self.frames["GALEX FUV"] * conversion_factor
        self.fuv_si.unit = "W/m2"

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

        ## GET THE GALAMETZ PARAMETERS
        a, b, c = self.get_galametz_parameters("MIPS 24mu", "Pacs blue", "Pacs red")

        assert a == 2.133
        assert b == 0.681
        assert c == 1.125

        # MIPS, PACS BLUE AND PACS RED CONVERTED TO LSUN (ABOVE)
        # Galametz (2013) formula for Lsun units
        tir_map = a * self.frames["MIPS 24mu"] + b * self.frames["Pacs blue"] + c * self.frames["Pacs red"]

        ## Convert the TIR map from Lsun to W / m2

        conversion_factor = 1.0

        # Conversion from Lsun to W

        conversion_factor *= solar_luminosity.to("W").value

        # Conversion from W [LUMINOSITY] to W / m2 [FLUX]
        distance = self.galaxy_properties.distance
        conversion_factor /= (4. * np.pi * distance**2).to("m2").value

        ## CONVERT AND SET NEW UNIT

        self.tir_si = Frame(tir_map * conversion_factor)
        self.tir_si.unit = "W/m2"

    # -----------------------------------------------------------------

    def get_galametz_parameters(self, *args):

        """
        This function ...
        :param args:
        :return:
        """

        galametz_column_names = {"MIPS 24mu": "c24",
                                 "Pacs blue": "c70",
                                 "Pacs green": "c100",
                                 "Pacs red": "c160",
                                 "SPIRE PSW": "c250"}

        # Needed column names
        needed_column_names = [galametz_column_names[filter_name] for filter_name in args]

        # List of not needed column names
        colnames = self.galametz.colnames
        not_needed_column_names = [name for name in colnames if name not in needed_column_names]

        not_needed_column_names.remove("R2")
        not_needed_column_names.remove("CV(RMSE)")

        # The parameters
        parameters = None

        # Loop over the entries in the galametz table
        for i in range(len(self.galametz)):

            if is_appropriate_galametz_entry(self.galametz, i, needed_column_names, not_needed_column_names):

                parameters = []
                for name in needed_column_names: parameters.append(self.galametz[name][i])
                break

            else: continue

        # Return the parameters
        return parameters

    # -----------------------------------------------------------------

    def make_tir_to_fuv(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the TIR to FUV map ...")

        # CALCULATE FUV AND TIR MAP IN W/M2 UNIT

        ## FUV IN W/M2

        ## TIR IN W/M2

        # CALCULATE TIR TO FUV RATIO

        # The ratio of TIR and FUV
        self.tir_to_fuv = self.tir_si / self.fuv_si
        self.log_tir_to_fuv = Frame(np.log10(self.tir_to_fuv))

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the TIR map
        self.write_tir()

        # Write the FUV map in SI units
        self.write_fuv()

        # Write the TIR to FUV ratio map
        self.write_tir_to_fuv()

        # Write the logarithm of TIR to FUV
        self.write_log_tir_to_fuv()

    # -----------------------------------------------------------------

    def write_tir(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the TIR map ...")

        # Determine path
        path = fs.join(self.maps_tirfuv_path, "TIR.fits")

        # Write
        self.tir_si.saveto(path)

    # -----------------------------------------------------------------

    def write_fuv(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the FUV map ...")

        # Determine path
        path = fs.join(self.maps_tirfuv_path, "FUV.fits")

        # Write
        self.fuv_si.saveto(path)

    # -----------------------------------------------------------------

    def write_tir_to_fuv(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the TIR to FUV map ...")

        # Determine path
        path = fs.join(self.maps_tirfuv_path, "TIR-FUV.fits")

        # Write
        self.tir_to_fuv.saveto(path)

    # -----------------------------------------------------------------

    def write_log_tir_to_fuv(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the logarithm of the TIR to FUV map ...")

        # Determine the path
        path = fs.join(self.maps_tirfuv_path, "logTIR-FUV.fits")

        # Write
        self.log_tir_to_fuv.saveto(path)

# -----------------------------------------------------------------

def is_appropriate_galametz_entry(table, i, needed_cols, not_needed_cols):

    """
    This function ...
    :return:
    """

    # Check if masked columns
    for name in not_needed_cols:

        # Verify that this entry is masked
        if not table[name].mask[i]: return False

    # Check needed cols
    for name in needed_cols:

        # Verify that this entry is not masked
        if table[name].mask[i]: return False

    # No mismatches found, thus appropriate entry
    return True

# -----------------------------------------------------------------
