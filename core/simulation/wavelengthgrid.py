#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.wavelengthgrid Contains the WavelengthGrid class.

# -----------------------------------------------------------------

# Import standard modules
import copy

# Import astronomical modules
from astropy.table import Table
from astropy.units import Unit

# Import the relevant PTS classes and modules
from ...core.tools import tables

# -----------------------------------------------------------------

class WavelengthGrid(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        """

        # Attributes
        #self.table = Table(names=["Wavelength", "Delta"], dtype=('f8', 'f8'))
        self.table = Table(names=["Wavelength"], dtype=('f8'))
        self.table["Wavelength"].unit = Unit("micron")
        #self.table["Delta"].unit = Unit("micron")

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Create a new class instance
        grid = cls()

        # Load the table
        grid.table = tables.from_file(path, format="ascii.ecsv")

        # Return the new instance
        return grid

    # -----------------------------------------------------------------

    @property
    def has_deltas(self):

        """
        This function ...
        :return:
        """

        return "Delta" in self.table.colnames

    # -----------------------------------------------------------------

    @classmethod
    def from_wavelengths(cls, wavelengths):

        """
        This function ...
        :param wavelengths:
        :return:
        """

        # Create a new class instance
        grid = cls()

        # Add the wavelengths
        grid.table["Wavelength"] = wavelengths

        # Return the new instance
        return grid

    # -----------------------------------------------------------------

    @classmethod
    def from_skirt_output(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Create a new class instance
        grid = cls()

        # Load the table
        table = tables.from_file(path, format="ascii")

        # Set the column names and units
        table.rename_column("col1", "Wavelength")
        table.rename_column("col2", "Delta")
        table["Wavelength"].unit = "micron"
        table["Delta"].unit = "micron"

        # Set the table
        grid.table = table

        # Return the new instance
        return grid

    # -----------------------------------------------------------------

    @classmethod
    def from_skirt_input(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        return NotImplementedError("Not implemented yet")

    # -----------------------------------------------------------------

    def __len__(self):

        """
        This function ...
        :return:
        """

        return len(self.table)

    # -----------------------------------------------------------------

    def add_point(self, wavelength):

        """
        This function ...
        :param wavelength:
        :return:
        """

        self.table.add_row([wavelength])

    # -----------------------------------------------------------------

    def wavelengths(self, unit=None, asarray=False, add_unit=True):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :return:
        """

        if asarray: return tables.column_as_array(self.table["Wavelength"], unit=unit)
        else: return tables.column_as_list(self.table["Wavelength"], unit=unit, add_unit=add_unit)

    # -----------------------------------------------------------------

    def deltas(self, unit=None, asarray=False, add_unit=True):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :return:
        """

        if asarray: return tables.column_as_array(self.table["Delta"], unit=unit)
        else: return tables.column_as_list(self.table["Delta"], unit=unit, add_unit=add_unit)

    # -----------------------------------------------------------------

    def to_skirt_input(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        table_copy = copy.deepcopy(self.table)

        # Write the wavelength grid
        table_copy.rename_column("Wavelength", str(len(self.table)))  # Trick to have the number of wavelengths in the first line (required for SKIRT)
        tables.write(table_copy, path, format="ascii")

    # -----------------------------------------------------------------

    def to_skirt_output(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        raise NotImplementedError("Not implemented yet")

    # -----------------------------------------------------------------

    def save(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Write the table
        tables.write(self.table, path, format="ascii.ecsv")

# -----------------------------------------------------------------
