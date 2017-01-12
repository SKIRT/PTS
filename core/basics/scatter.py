#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.scatter Contains the Scatter class and derived classes.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from .table import SmartTable
from ..tools.strings import alphabet
from .range import RealRange
from ..tools import tables

# -----------------------------------------------------------------

class Scatter(SmartTable):

    """
    This class ...
    """

    column_info = []

    # -----------------------------------------------------------------

    @classmethod
    def initialize(cls, variables, units, descriptions):

        """
        This function ...
        :param variables: list of the names of the variables   [OR: AN INTER WITH THE NUMBER OF VARIABLES]
        :param units: a dict with units  [OR: A LIST WITH A LENGTH OF THE NUMBER OF VARIABLES]
        :param descriptions: a dict with descriptions [OR: A LIST WITH A LENGTH OF THE NUMBER OF VARIABLES]
        :return:
        """

        if isinstance(variables, int):

            # Loop over the variables
            for index in range(variables):

                # Add column
                cls.column_info.append((alphabet[index], float, units[index], descriptions[index]))

        else:

            # Loop over the variables
            for variable in variables:

                unit = units[variable] if variable in units else None
                description = descriptions[variable] if variable in descriptions else None

                # Add column
                cls.column_info.append((variable, float, unit, description))

        # Call the initialize function of the SmartTable table function
        return super(Scatter, cls).initialize()

    # -----------------------------------------------------------------

    @property
    def nvariables(self):

        """
        This function ...
        :return:
        """

        return len(self.colnames)

    # -----------------------------------------------------------------

    def add_point(self, *args):

        """
        This function ...
        :param args: the values
        :return:
        """

        # Check length
        if len(args) != self.nvariables: raise ValueError("The number of variables is " + str(self.nvariables))

        # Add a row
        self.add_row(args)

    # -----------------------------------------------------------------

    def range_of(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        values = tables.column_as_array(self[name], unit=self[name].unit)
        mask = np.logical_not(np.isnan(values))
        values = values[mask]

        range = RealRange(np.min(values), np.max(values))
        return range

    # -----------------------------------------------------------------

    def ranges(self, as_list=False):

        """
        This function ...
        :param as_list:
        :return:
        """

        result = None

        if as_list:

            result = []
            for name in self.colnames: result.append(self.range_of(name))

        else:

            result = dict()
            for name in self.colnames: result[name] = self.range_of(name)

        # Return
        return result

# -----------------------------------------------------------------

class Scatter3D(Scatter):

    """
    This class ...
    """

    @classmethod
    def initialize(cls, x_label, y_label, z_label, x_unit, y_unit, z_unit, x_description, y_description, z_description):

        """
        This function ...
        :param x_label:
        :param y_label:
        :param z_label:
        :param x_unit:
        :param y_unit:
        :param z_unit:
        :param x_description:
        :param y_description:
        :param z_description:
        :return:
        """

        # Call the initialize function of the base class
        return super(Scatter3D, cls).initialize([x_label, y_label, z_label], {x_label: x_unit, y_label: y_unit, z_label: z_unit}, {x_label: x_description, y_label: y_description, z_label: z_description})

    # -----------------------------------------------------------------

    def add_point(self, x, y, z):

        """
        This function ...
        :param x:
        :param y:
        :param z:
        :return:
        """

        super(Scatter3D, self).add_point(x, y ,z)

    # -----------------------------------------------------------------

    @property
    def x_label(self):

        """
        This function ...
        :return:
        """

        return self.colnames[0]

    # -----------------------------------------------------------------

    @property
    def y_label(self):

        """
        This function ...
        :return:
        """

        return self.colnames[1]

    # -----------------------------------------------------------------

    @property
    def z_label(self):

        """
        This fucntion ...
        :return:
        """

        return self.colnames[2]

    # -----------------------------------------------------------------

    @property
    def x_range(self):

        """
        This function ...
        :return:
        """

        return self.range_of(self.x_label)

    # -----------------------------------------------------------------

    @property
    def y_range(self):

        """
        This function ...
        :return:
        """

        return self.range_of(self.y_label)

    # -----------------------------------------------------------------

    @property
    def z_range(self):

        """
        This function ...
        :return:
        """

        return self.range_of(self.z_label)

    # -----------------------------------------------------------------

    def x(self, unit=None, asarray=False, add_unit=True):

        """
        This function ...
        :return:
        """

        if asarray: return tables.column_as_array(self[self.x_label], unit=unit)
        else: return tables.column_as_list(self[self.x_label], unit=unit, add_unit=add_unit)

    # -----------------------------------------------------------------

    def y(self, unit=None, asarray=False, add_unit=True):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :return:
        """

        if asarray: return tables.column_as_array(self[self.y_label], unit=unit)
        else: return tables.column_as_list(self[self.y_label], unit=unit, add_unit=add_unit)

    # -----------------------------------------------------------------

    def z(self, unit=None, asarray=False, add_unit=True):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :return:
        """

        if asarray: return tables.column_as_array(self[self.z_label], unit=unit)
        else: return tables.column_as_list(self[self.z_label], unit=unit, add_unit=add_unit)

# -----------------------------------------------------------------
