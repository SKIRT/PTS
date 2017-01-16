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

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(FluxErrorTable, self).__init__(*args, **kwargs)

        # Add column info
        self.add_column_info("Instrument", str, None, "name of the instrument")
        self.add_column_info("Band", str, None, "name of the band")
        self.add_column_info("Calibration error-", float, None, "calibration error lower")
        self.add_column_info("Calibration error+", float, None, "calibration error upper")
        self.add_column_info("Aperture noise-", float, None, "aperture noise lower")
        self.add_column_info("Aperture noise+", float, None, "aperture noise upper")
        self.add_column_info("Total error-", float, None, "total error lower")
        self.add_column_info("Total error+", float, None, "total error upper")
        self.add_column_info("Total relative error", float, None, "total relative error")

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

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Check
        if "labels" in kwargs: from_astropy = False
        else: from_astropy = True

        # Get properties
        if not from_astropy: labels = kwargs.pop("labels")
        else: labels = None

        # Call the constructor of the base class
        super(FluxDifferencesTable, self).__init__(*args, **kwargs)

        # Add column info
        if not from_astropy:

            # Add column info
            self.add_column_info("Instrument", str, None, "name of of the instrument")
            self.add_column_info("Band", str, None, "name of the band")
            self.add_column_info("Flux", float, None, "the flux value")

            # Loop over the labels
            for label in labels: self.add_column_info(label + " relative difference", float, None, "relative difference with the " + label + " flux")

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
