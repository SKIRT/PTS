#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.photometry.tables Contains the FluxErrorTable class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.table import SmartTable

# -----------------------------------------------------------------

class FluxErrorTable(SmartTable):
    
    """
    This class...
    """

    column_info = [("Instrument", str, None, "name of the instrument"),
                   ("Band", str, None, "name of the band"),
                   ("Calibration error-", float, None, "calibration error lower"),
                   ("Calibration error+", float, None, "calibration error upper"),
                   ("Aperture noise-", float, None, "aperture noise lower"),
                   ("Aperture noise+", float, None, "aperture noise upper"),
                   ("Total error-", float, None, "total error lower"),
                   ("Total error+", float, None, "total error upper")]

    # -----------------------------------------------------------------

    def add_entry(self, fltr, calibration_error, aperture_noise, total_error):

        """
        This function ...
        :param fltr:
        :param calibration_error:
        :param aperture_noise:
        :param total_error:
        :return:
        """

        # Set the values
        values = [fltr.instrument, fltr.band, calibration_error.lower, calibration_error.upper, aperture_noise.lower, aperture_noise.upper, total_error.lower, total_error.upper]

        # Add a row to the table
        self.add_row(values)

# -----------------------------------------------------------------
