#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.sed Contains the SED class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy import units as u
from astropy.table import Table

# Import the relevant PTS classes and modules
from ...core.tools import tables
from ...core.basics.errorbar import ErrorBar

# -----------------------------------------------------------------

class ObservedSED(object):

    """
    This function ...
    """

    def __init__(self):

        """
        This function ...
        :return:
        """

        # Attributes
        self.table = Table(names=["Observatory", "Instrument", "Band", "Wavelength", "Flux", "Error-", "Error+"], dtype=('S10', 'S10', 'S10', 'f8', 'f8', 'f8', 'f8'))
        self.table["Wavelength"].unit = u.Unit("micron")
        self.table["Flux"].unit = u.Unit("Jy")
        self.table["Error-"].unit = u.Unit("Jy")
        self.table["Error+"].unit = u.Unit("Jy")

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Create a new observed SED
        sed = cls()

        # Open the SED table
        sed.table = tables.from_file(path, format="ascii.commented_header")
        sed.table["Wavelength"].unit = u.Unit("micron")
        sed.table["Flux"].unit = u.Unit("Jy")
        sed.table["Error-"].unit = u.Unit("Jy")
        sed.table["Error+"].unit = u.Unit("Jy")

        # Return the observed SED
        return sed

    # -----------------------------------------------------------------

    def add_entry(self, filter, flux, error):

        """
        This function ...
        :param filter:
        :param flux:
        :param error:
        :return:
        """

        self.table.add_row([filter.observatory, filter.instrument, filter.band, filter.pivotwavelength(), flux, error.lower, error.upper])

    # -----------------------------------------------------------------

    def instruments(self):

        """
        This function ...
        :return:
        """

        return tables.column_as_list(self.table["Instrument"])

    # -----------------------------------------------------------------

    def bands(self):

        """
        This function ...
        :return:
        """

        return tables.column_as_list(self.table["Band"])

    # -----------------------------------------------------------------

    def wavelengths(self, unit=None):

        """
        This function ...
        :param unit:
        :return:
        """

        return tables.column_as_list(self.table["Wavelength"], unit=unit)

    # -----------------------------------------------------------------

    def fluxes(self, unit=None):

        """
        This function ...
        :param unit:
        :return:
        """

        return tables.column_as_list(self.table["Flux"], unit=unit)

    # -----------------------------------------------------------------

    def errors(self, unit=None):

        """
        This function ...
        :param unit:
        :return:
        """

        return tables.columns_as_objects([self.table["Error-"], self.table["Error+"]], ErrorBar, unit=unit)

    # -----------------------------------------------------------------

    def save(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Sort the table by wavelength
        self.table.sort("Wavelength")

        # Write the observed SED
        tables.write(self.table, path)

# -----------------------------------------------------------------

class SED(object):
    
    """
    This class...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Attributes
        self.table = None

    # -----------------------------------------------------------------

    def wavelengths(self, unit=None):

        """
        This function ...
        :param unit:
        :return:
        """

        return tables.column_as_list(self.table["Wavelength"], unit=unit)

    # -----------------------------------------------------------------

    def fluxes(self, unit=None):

        """
        This function ...
        :param unit:
        :return:
        """

        return tables.column_as_list(self.table["Flux"], unit=unit)

    # -----------------------------------------------------------------

    def errors(self, unit=None):

        """
        This function ...
        :param unit:
        :return:
        """

        return tables.columns_as_objects([self.table["Error-"], self.table["Error+"]], ErrorBar, unit=unit)

    # -----------------------------------------------------------------

    @property
    def has_errors(self):

        """
        This function ...
        :return:
        """

        return "Error-" in self.table.colnames and "Error+" in self.table.colnames

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Create a new SED
        sed = cls()

        # From SKIRT:
        # column 1: lambda (micron)
        # column 2: total flux; lambda*F_lambda (W/m2)

        from ..preparation import unitconversion

        # Open the SED table
        sed.table = tables.from_file(path, format="ascii.no_header")
        sed.table.rename_column("col1", "Wavelength")
        sed.table.rename_column("col2", "Flux")
        sed.table["Wavelength"].unit = u.Unit("micron")

        jansky_column = []

        for i in range(len(sed.table)):

            # Get the flux density in W / m2 and the wavelength in micron
            neutral_fluxdensity = sed.table["Flux"][i]
            wavelength = sed.table["Wavelength"][i] * u.Unit("micron")

            # Convert to Jansky
            jansky = unitconversion.neutral_fluxdensity_to_jansky(neutral_fluxdensity, wavelength)

            # Add the fluxdensity in Jansky to the new column
            jansky_column.append(jansky)

        # Add the flux column in Jansky
        sed.table.remove_column("Flux")
        sed.table["Flux"] = jansky_column
        sed.table["Flux"].unit = "Jy"

        # Return the SED
        return sed

    # -----------------------------------------------------------------

    def save(self, path):

        """
        This function ...
        :param path
        :return:
        """

        # Write the SED table to file
        tables.write(self.table, path)

# -----------------------------------------------------------------
