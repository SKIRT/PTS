#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.heating.tables Contains different table classes.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....core.basics.table import SmartTable
from ....core.tools import arrays
from ....core.units.unit import PhotometricUnit

# -----------------------------------------------------------------

class AbsorptionTable(SmartTable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(AbsorptionTable, self).__init__(*args, **kwargs)

        # Add columns
        #self.add_column_info("Run name", str, None, "Name of the analysis run")
        #self.add_column_info("Host id", str, None, "Cache remote host ID")

        # X coordinate of cell center
        # Y coordinate of cell center
        # Z coordinate of cell center

        # Absorbed bolometric luminosity of the total stellar population
        # Absorbed bolometric luminosity of the old stellar population
        # Absorbed bolometric luminosity of the young stellar population
        # Absorbed bolometric luminosity of the ionizing stellar population

        # Add column info
        self.add_column_info("x", float, "pc", "X coordinate of cell center")
        self.add_column_info("y", float, "pc", "Y coordinate of cell center")
        self.add_column_info("z", float, "pc", "Z coordinate of cell center")
        self.add_column_info("total", float, "W", "absorbed bolometric luminosity of the total stellar population")
        self.add_column_info("old", float, "W", "absorbed bolometric luminosity of the old stellar population")
        self.add_column_info("young", float, "W", "absorbed bolometric luminosity of the young stellar population")
        self.add_column_info("ionizing", float, "W", "absorbed bolometric luminosity of the ionizing stellar population")

    # -----------------------------------------------------------------

    def add_entry(self, x, y, z, total, old, young, ionizing):

        """
        This function ...
        :param x:
        :param y:
        :param z:
        :param total:
        :param old:
        :param young:
        :param ionizing:
        :return:
        """

        values = [x, y, z, total, old, young, ionizing]
        self.add_row(values)

    # -----------------------------------------------------------------

    def total(self, unit="W", add_unit=True, asarray=False):

        """
        This function ...
        :param unit:
        :param add_unit:
        :param asarray:
        :return:
        """

        if asarray: return arrays.plain_array(self["total"], unit=unit, array_unit=self["total"].unit)
        else: return arrays.array_as_list(self["total"], unit=unit, add_unit=add_unit, array_unit=self["total"].unit)

    # -----------------------------------------------------------------

    def old(self, unit="W", add_unit=True, asarray=False):

        """
        This function ...
        :param unit:
        :param add_unit:
        :param asarray:
        :return:
        """

        if asarray: return arrays.plain_array(self["old"], unit=unit, array_unit=self["old"].unit)
        else: return arrays.array_as_list(self["old"], unit=unit, add_unit=add_unit, array_unit=self["old"].unit)

    # -----------------------------------------------------------------

    def young(self, unit="W", add_unit=True, asarray=False):

        """
        This function ...
        :param unit:
        :param add_unit:
        :param asarray:
        :return:
        """

        if asarray: return arrays.plain_array(self["young"], unit=unit, array_unit=self["young"].unit)
        else: return arrays.array_as_list(self["young"], unit=unit, add_unit=add_unit, array_unit=self["young"].unit)

    # -----------------------------------------------------------------

    def ionizing(self, unit="W", add_unit=True, asarray=False):

        """
        Thisn function ...
        :param unit:
        :param add_unit:
        :param asarray:
        :return:
        """

        if asarray: return arrays.plain_array(self["ionizing"], unit=unit, array_unit=self["ionizing"].unit)
        else: return arrays.array_as_list(self["ionizing"], unit=unit, add_unit=add_unit, array_unit=self["ionizing"].unit)

    # -----------------------------------------------------------------

    def unevolved(self, unit="W", add_unit=True, asarray=True):

        """
        This function ...
        :param unit:
        :param add_unit:
        :param asarray:
        :return:
        """

        # Create data
        young = self.young(asarray=True, unit=unit)
        ionizing = self.ionizing(asarray=True, unit=unit)
        unevolved = young + ionizing

        # Create unit object
        unit = PhotometricUnit(unit)

        # Return
        if asarray: return unevolved
        elif add_unit: [value * unit for value in unevolved]
        else: return list(unevolved)

# -----------------------------------------------------------------
