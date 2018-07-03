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
from ..tools import tables, arrays
from .relation import Relation

# -----------------------------------------------------------------

class Scatter2D(Relation):

    """
    This class ...
    """

# -----------------------------------------------------------------

class Scatter(SmartTable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Check
        if "variables" in kwargs: from_astropy = False
        else: from_astropy = True

        # Get properties
        if not from_astropy:
            variables = kwargs.pop("variables")
            units = kwargs.pop("units")
            descriptions = kwargs.pop("descriptions")
        else: variables = units = descriptions = None

        # Call the constructor of the base class
        super(Scatter, self).__init__(*args, **kwargs)

        # Add column info
        if not from_astropy:

            if isinstance(variables, int):

                # Loop over the variables
                for index in range(variables): self.add_column_info(alphabet[index], float, units[index], descriptions[index])

            else:

                # Loop over the variables
                for variable in variables:

                    unit = units[variable] if variable in units else None
                    description = descriptions[variable] if variable in descriptions else None

                    # Add column
                    self.add_column_info(variable, float, unit, description)

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

        values = arrays.plain_array(self[name], unit=self[name].unit)
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

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        if "x_label" in kwargs: from_astropy = False
        else: from_astropy = True

        # Get properties
        if not from_astropy:

            x_label = kwargs.pop("x_label")
            y_label = kwargs.pop("y_label")
            z_label = kwargs.pop("z_label")
            x_unit = kwargs.pop("x_unit")
            y_unit = kwargs.pop("y_unit")
            z_unit = kwargs.pop("z_unit")
            x_description = kwargs.pop("x_description")
            y_description = kwargs.pop("y_description")
            z_description = kwargs.pop("z_description")

            kwargs["variables"] = [x_label, y_label, z_label]
            kwargs["units"] = {x_label: x_unit, y_label: y_unit, z_label: z_unit}
            kwargs["descriptions"] = {x_label: x_description, y_label: y_description, z_label: z_description}

        else: x_label = y_label = z_label = x_unit = y_unit = z_unit = x_description = y_description = z_description = None

        # Call the constructor of the base class
        super(Scatter3D, self).__init__(*args, **kwargs)

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

        if asarray: return arrays.plain_array(self[self.x_label], unit=unit)
        else: return arrays.array_as_list(self[self.x_label], unit=unit, add_unit=add_unit)

    # -----------------------------------------------------------------

    def y(self, unit=None, asarray=False, add_unit=True):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :return:
        """

        if asarray: return arrays.plain_array(self[self.y_label], unit=unit)
        else: return arrays.array_as_list(self[self.y_label], unit=unit, add_unit=add_unit)

    # -----------------------------------------------------------------

    def z(self, unit=None, asarray=False, add_unit=True):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :return:
        """

        if asarray: return arrays.plain_array(self[self.z_label], unit=unit)
        else: return arrays.array_as_list(self[self.z_label], unit=unit, add_unit=add_unit)

# -----------------------------------------------------------------
