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
                   ("Total error+", float, None, "total error upper"),
                   ("Total relative error", float, None, "total relative error")]

    # -----------------------------------------------------------------

    def add_entry(self, fltr, calibration_error, aperture_noise, total_error, total_relative_error):

        """
        This function ...
        :param fltr:
        :param calibration_error:
        :param aperture_noise:
        :param total_error:
        :param total_relative_error:
        :return:
        """

        # Set the values
        values = [fltr.instrument, fltr.band, calibration_error.lower, calibration_error.upper, aperture_noise.lower, aperture_noise.upper, total_error.lower, total_error.upper, total_relative_error]

        # Add a row to the table
        self.add_row(values)

# -----------------------------------------------------------------

class FluxDifferencesTable(SmartTable):

    """
    This class ...
    """

    column_info = [("Instrument", str, None, "name of of the instrument"),
                   ("Band", str, None, "name of the band"),
                   ("Flux", float, None, "the flux value")]

    # -----------------------------------------------------------------

    @classmethod
    def initialize(cls, labels):

        """
        This function ...
        :param labels:
        :return:
        """

        # Add the column info
        for label in labels: cls.column_info.append((label + " relative difference", float, None, "relative difference with the " + label + " flux"))

        # Call the initialize function of the generations table function
        return super(FluxDifferencesTable, cls).initialize()

    # -----------------------------------------------------------------

    @property
    def reference_labels(self):

        """
        This function ...
        :return:
        """

        labels = []

        for name in self.colnames:

            label = name.split(" relative difference")[0]
            labels.append(label)

        # Return the list of labels
        return labels

    # -----------------------------------------------------------------

    def add_entry(self, fltr, flux, flux_differences):

        """
        This function ...
        :param fltr:
        :param flux:
        :param flux_differences:
        :return:
        """

        values = [fltr.instrument, fltr.band, flux]

        # Add the parameter values
        for label in self.reference_labels: values.append(flux_differences[label])

        # Add a row to the table
        self.add_row(values)

# -----------------------------------------------------------------
