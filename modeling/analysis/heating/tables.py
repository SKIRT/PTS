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

# Import standard modules
from collections import OrderedDict

# Import astronomical modules
from astropy.table import Column

# Import the relevant PTS classes and modules
from ....core.basics.table import SmartTable
from ....core.tools import arrays
from ....core.units.unit import PhotometricUnit
from ....magic.basics.coordinate import PhysicalCoordinate, PhysicalCoordinate3D
from ....magic.basics.vector import Position, Position3D

# -----------------------------------------------------------------

class AbsorptionTable(SmartTable):

    """
    This class ...
    """

    # Add column info
    _column_info = OrderedDict()
    _column_info["x"] = (float, "pc", "X coordinate of cell center")
    _column_info["y"] = (float, "pc", "Y coordinate of cell center")
    _column_info["z"] = (float, "pc", "Z coordinate of cell center")
    _column_info["total"] = (float, "W", "absorbed bolometric luminosity of the total stellar population")
    _column_info["old"] = (float, "W", "absorbed bolometric luminosity of the old stellar population")
    _column_info["young"] = (float, "W", "absorbed bolometric luminosity of the young stellar population")
    _column_info["ionizing"] = (float, "W", "absorbed bolometric luminosity of the ionizing stellar population")
    _column_info["extra"] = (float, "W", "absorbed bolometric luminosity of the extra component")

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(AbsorptionTable, self).__init__(*args, **kwargs)

        # Add column info
        self.add_all_column_info(self._column_info)

    # -----------------------------------------------------------------

    @classmethod
    def from_columns(cls, x, y, z, total, old, young, ionizing, extra):

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

        new = cls()
        new._setup()

        new.remove_columns(["x", "y", "z", "total", "old", "young", "ionizing", "extra"])

        if not isinstance(x, Column): x = Column(data=x)
        if not isinstance(y, Column): y = Column(data=y)
        if not isinstance(z, Column): z = Column(data=z)
        if not isinstance(total, Column): total = Column(data=total)
        if not isinstance(old, Column): old = Column(data=old)
        if not isinstance(young, Column): young = Column(data=young)
        if not isinstance(ionizing, Column): ionizing = Column(data=ionizing)
        if not isinstance(extra, Column): extra = Column(data=extra)

        new.add_columns([x, y, z], copy=False, names=["x", "y", "z"])
        #new.rename_column("X coordinate of cell center", "x")
        #new.rename_column("Y coordinate of cell center", "y")
        #new.rename_column("Z coordinate of cell center", "z")
        #new.rename_column(x.name, "x")
        #new.rename_column(y.name, "y")
        #new.rename_column(z.name, "z")
        new["x"].unit = "pc"
        new["y"].unit = "pc"
        new["z"].unit = "pc"

        new.add_column(total, name="total")
        #new.rename_column("Absorbed bolometric luminosity", "total")
        #new.rename_column(total.name, "total")
        new["total"].unit = "W"

        new.add_column(old, name="old")
        #new.rename_column("Absorbed bolometric luminosity", "old")
        #new.rename_column(old.name, "old")
        new["old"].unit = "W"

        new.add_column(young, name="young")
        #new.rename_column("Absorbed bolometric luminosity", "young")
        #new.rename_column(young.name, "young")
        new["young"].unit = "W"

        new.add_column(ionizing, name="ionizing")
        #new.rename_column("Absorbed bolometric luminosity", "ionizing")
        #new.rename_column(ionizing.name, "ionizing")
        new["ionizing"].unit = "W"

        new.add_column(ionizing, name="extra")
        # new.rename_column("Absorbed bolometric luminosity", "extra")
        # new.rename_column(extra.name, "extra")
        new["extra"].unit = "W"

        # Return the new table
        return new

    # -----------------------------------------------------------------

    def add_entry(self, x, y, z, total, old, young, ionizing, extra):

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

        values = [x, y, z, total, old, young, ionizing, extra]
        self.add_row(values)

    # -----------------------------------------------------------------

    @property
    def ncells(self):

        """
        This function ...
        :return:
        """

        return len(self)

    # -----------------------------------------------------------------

    def planar_coordinates(self, unit="pc", add_unit=True):

        """
        Thisf unction ...
        :param unit:
        :param add_unit:
        :return:
        """

        coordinates = []

        # Loop over the cells
        for index in range(self.ncells):

            # Get x and y
            x = self["x"][index] * self["x"].unit
            y = self["y"][index] * self["y"].unit

            # Convert
            x = x.to(unit)
            y = y.to(unit)

            # With unit
            if add_unit: coordinate = PhysicalCoordinate(x, y)

            # Without unit
            else: coordinate = Position(x, y)

            # Add the coordinate
            coordinates.append(coordinate)

        # Return
        return coordinates

    # -----------------------------------------------------------------

    def coordinates(self, unit="pc", add_unit=True):

        """
        This function ...
        :param unit:
        :param add_unit:
        :return:
        """

        coordinates = []

        # Loop over the cells
        for index in range(self.ncells):

            # Get x and y
            x = self["x"][index] * self["x"].unit
            y = self["y"][index] * self["y"].unit
            z = self["z"][index] * self["z"].unit

            # Convert
            x = x.to(unit)
            y = y.to(unit)
            z = z.to(unit)

            # With unit
            if add_unit: coordinate = PhysicalCoordinate3D(x, y, z)

            # Without unit
            else: coordinate = Position3D(x, y, z)

            # Add the coordinate
            coordinates.append(coordinate)

        # Return
        return coordinates

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

    def extra(self, unit="W", add_unit=True, asarray=False):

        """
        This function ...
        :param unit:
        :param add_unit:
        :param asarray:
        :return:
        """

        if asarray: return arrays.plain_array(self["extra"], unit=unit, array_unit=self["extra"].unit)
        else: return arrays.array_as_list(self["extra"], unit=unit, add_unit=add_unit, array_unit=self["extra"].unit)

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
