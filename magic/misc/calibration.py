#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.misc.calibration Contains the CalibrationError class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# -----------------------------------------------------------------

calibration_errors = {"GALEX FUV": "0.05 mag",
                      "GALEX NUV": "0.03 mag",
                      "UVOT UVW2": "2.18%",
                      "UVOT UVM2": "0.63%",
                      "UVOT UVW1": "3.09%",
                      "Mosaic Halpha": "5%",
                      "Halpha": "5%",
                      "SDSS u": "1.3%", # source: http://www.sdss.org/dr12/scope/
                      "SDSS g": "0.8%", # source: http://www.sdss.org/dr12/scope/
                      "SDSS r": "0.8%", # source: http://www.sdss.org/dr12/scope/
                      "SDSS i": "0.7%", # source: http://www.sdss.org/dr12/scope/
                      "SDSS z": "0.8%", # source: http://www.sdss.org/dr12/scope/
                      "2MASS J": "0.03 mag",
                      "2MASS H": "0.03 mag",
                      "2MASS Ks": "0.03 mag",
                      "WISE W1": "2.4%",
                      "IRAC I1": "1.8%",
                      "IRAC I2": "1.9%",
                      "WISE W2": "2.8%",
                      "IRAC I3": "2.0%",
                      "IRAC I4": "2.1%",
                      "WISE W3": "4.5%",
                      "WISE W4": "5.7%",
                      "MIPS 24mu": "4%",
                      "MIPS 70mu": "5%",
                      # ref: Absolute_Calibration_and_Characterization_of_the_Multiband_Imaging_Photometer_for_Spitzer._II._70_micron_Imaging
                      "MIPS 160mu": "5%",
                      "Pacs blue": "7%",
                      "Pacs green:": "7%",
                      "Pacs red": "7%",
                      "SPIRE PSW": "4%",
                      "SPIRE PMW": "4%",
                      "SPIRE PLW": "4%"}

# -----------------------------------------------------------------

class CalibrationError(object):

    """
    This class ...
    """

    def __init__(self, value, unit):

        """
        The constructor ...
        :param value:
        :param unit:
        """

        # value and unit
        self.value = value
        self.unit = unit

    # -----------------------------------------------------------------

    @classmethod
    def from_string(cls, string):

        """
        This function ...
        :param string:
        :return:
        """

        if "mag" in string:
            value = float(string.split(" mag")[0])
            unit = "mag"
        elif "%" in string:
            value = float(string.split("%")[0])
            unit = "%"
        else: raise ValueError("An error occured")

        # Return a new class instance
        return cls(value, unit)

    # -----------------------------------------------------------------

    @classmethod
    def from_filter_name(cls, filter_name):

        """
        This function ...
        :param filter_name:
        :return:
        """

        # Set the calibration error
        calibration_error = calibration_errors[filter_name]

        # Create and return the calibration error instance
        return cls.from_string(calibration_error)

    # -----------------------------------------------------------------

    @classmethod
    def from_filter(cls, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        # Get the filter name
        filter_name = str(fltr)

        # Return the calibration error
        return cls.from_filter_name(filter_name)

    # -----------------------------------------------------------------

    @classmethod
    def from_instrument_and_band(cls, instrument, band):

        """
        This function ...
        :param instrument:
        :param band:
        :return:
        """

        # Get the filter name
        filter_name = instrument + " " + band

        # Return the calibration error
        return cls.from_filter_name(filter_name)

    # -----------------------------------------------------------------

    @property
    def magnitude(self):

        """
        This function ...
        :return:
        """

        if self.unit == "mag": return True
        else: return False

    # -----------------------------------------------------------------

    @property
    def percentage(self):

        """
        This function ...
        :return:
        """

        if self.unit == "%": return True
        else: return False

# -----------------------------------------------------------------

