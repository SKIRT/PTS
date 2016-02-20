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
        self.table = Table(names=["Observatory", "Instrument", "Band", "Wavelength", "Flux", "Error"], dtype=('S10', 'S10', 'S10', 'f8', 'f8', 'f8'))
        self.table["Wavelength"].unit = u.Unit("micron")
        self.table["Flux"].unit = u.Unit("Jy")

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls):

        """
        This function ...
        :return:
        """

        # Create a new observed SED
        sed = cls()

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

        self.table.add_row([filter.observatory, filter.instrument, filter.band, filter.pivotwavelength(), flux, error])

    # -----------------------------------------------------------------

    @property
    def instruments(self):

        """
        This function ...
        :return:
        """

        return self.table["Instrument"]

    # -----------------------------------------------------------------

    @property
    def bands(self):

        """
        This function ...
        :return:
        """

        return self.table["Band"]

    # -----------------------------------------------------------------

    @property
    def wavelengths(self):

        """
        This function ...
        :return:
        """

        return self.table["Wavelength"]

    # -----------------------------------------------------------------

    @property
    def fluxes(self):

        """
        This function ...
        :return:
        """

        return self.table["Flux"]

    # -----------------------------------------------------------------

    @property
    def errors(self):

        """
        This function ...
        :return:
        """

        return self.table["Error"]

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
        self.errors = []

    # -----------------------------------------------------------------

    @property
    def wavelengths(self):

        """
        This function ...
        :return:
        """

        return self.table["Wavelength"]

    # -----------------------------------------------------------------

    @property
    def fluxes(self):

        """
        This function ...
        :return:
        """

        return self.table["Flux"]

    # -----------------------------------------------------------------

    @property
    def has_errors(self):

        """
        This function ...
        :return:
        """

        return len(self.errors) > 0 and len(self.errors) == len(self.table["Wavelength"])

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

        # Open the SED table
        sed.table = tables.from_file(path, format="ascii.no_header")
        sed.table.rename_column("col1", "Wavelength")
        sed.table.rename_column("col2", "Flux")
        sed.table["Wavelength"].unit = u.Unit("micron")
        sed.table["Flux"].unit = u.Unit("W/m2")

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
