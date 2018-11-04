#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.relation Contains the Relation class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
import copy

# Import astronomical modules
from astropy.table import MaskedColumn

# Import the relevant PTS classes and modules
from .table import SmartTable
from ..tools import arrays

# -----------------------------------------------------------------

errormin_name = "Error-"
errorplus_name = "Error+"

# -----------------------------------------------------------------

class Relation(SmartTable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        if kwargs.get("from_astropy", None) is None:
            if "x_unit" in kwargs: from_astropy = False
            else: from_astropy = True
        else: from_astropy = kwargs.pop("from_astropy")

        # Get properties
        if not from_astropy:

            x_unit = kwargs.pop("x_unit", None)
            y_unit = kwargs.pop("y_unit", None)

            x_name = kwargs.pop("x_name", "x")
            if x_name is None: x_name = "x"

            y_name = kwargs.pop("y_name", "y")
            if y_name is None: y_name = "y"

            x_description = kwargs.pop("x_description", "x values")
            if x_description is None: x_description = "x values"

            y_description = kwargs.pop("y_description", "y values")
            if y_description is None: y_description = "y values"

            names = kwargs.pop("names", None)
            units = kwargs.pop("units", None)

        else: x_unit = y_unit = x_name = y_name = x_description = y_description = names = units = None

        # Call the constructor of the base class
        super(Relation, self).__init__(*args, **kwargs)

        # Set properties
        if not from_astropy:

            # Set the column info
            self.add_column_info(x_name, float, x_unit, x_description)
            self.add_column_info(y_name, float, y_unit, y_description)

            # Set x name and y name
            self.x_name = x_name
            self.y_name = y_name

            # Has auxilary axes?
            #names = kwargs.pop("names", None)
            #units = kwargs.pop("units", None)
            #print(names, units)

            # Auxilary axes?
            if names is not None and len(names) > 2:

                aux_names = names[2:]
                naux = len(aux_names)
                #print(units)
                if units is not None:
                    from . import containers
                    if isinstance(units, dict): aux_units = containers.sequence_from_dict(units, aux_names, none_if_not_contains=True)
                    else: aux_units = units[2:]
                else: aux_units = [None] * naux
                #aux_units = units[2:] if units is not None else [None] * naux

                # Add column info for auxilary axes
                for name, unit in zip(aux_names, aux_units): self.add_column_info(name, float, unit, "auxilary axis")

    # -----------------------------------------------------------------

    def copy(self, aux=True):

        """
        This function ...
        :param aux:
        :return:
        """

        new = copy.deepcopy(self)
        if hasattr(self, "x_name"): new.x_name = self.x_name
        if hasattr(self, "y_name"): new.y_name = self.y_name
        return new

    # -----------------------------------------------------------------

    @classmethod
    def from_columns(cls, *columns, **kwargs):

        """
        This function ...
        :param columns:
        :param kwargs:
        :return:
        """

        kwargs["from_astropy"] = False

        # x and y name
        if kwargs.get("names", None) is not None:

            names = kwargs.get("names")
            x_name = names[0]
            y_name = names[1]

            # Set
            kwargs["x_name"] = x_name
            kwargs["y_name"] = y_name

        #else: x_name, y_name = "x", "y"

        # x and y unit
        if kwargs.get("units", None) is not None:
            units = kwargs.get("units")

            from ..tools import types
            if types.is_sequence_or_tuple(units):
                x_unit = units[0]
                y_unit = units[1]
            elif types.is_dictionary(units):
                if "x_name" not in kwargs: raise ValueError("Cannot determine x unit if x name is not explicitely defined")
                if "y_name" not in kwargs: raise ValueError("Cannot determine y unit if y name is not explicitely defined")
                x_unit = units[kwargs["x_name"]] if kwargs["x_name"] in units else None
                y_unit = units[kwargs["y_name"]] if kwargs["y_name"] in units else None
            else: raise ValueError("Invalid type for 'units': must be sequence or dictionary")

            # Set
            kwargs["x_unit"] = x_unit
            kwargs["y_unit"] = y_unit

        #else: x_unit = y_unit = None

        # Add auxilary axes?
        aux = kwargs.pop("aux", None)
        if aux is not None:
            from .containers import sequence_from_dict
            naux = len(aux)
            columns = list(columns) + aux.values()  # columns = [x, y] + aux.values()
            column_names = [kwargs["x_name"], kwargs["y_name"]] + aux.keys()
            aux_units = kwargs.pop("aux_units", None)
            if aux_units is None: aux_units = [None] * naux
            else: aux_units = sequence_from_dict(aux_units, aux.keys(), none_if_not_contains=True) # in the order of the aux columns
            column_units = [kwargs["x_unit"], kwargs["y_unit"]] + aux_units
            kwargs["names"] = column_names
            kwargs["units"] = column_units
        # else: #columns = [x, y]
        if "aux_units" in kwargs: kwargs.pop("aux_units")

        # x and y descriptions
        #if kwargs.get("descriptions", None) is not None:
        descriptions = kwargs.get("descriptions", ("no description", "no description"))
        x_description = descriptions[0]
        y_description = descriptions[1]
        kwargs["x_description"] = x_description
        kwargs["y_description"] = y_description

        #print(kwargs)

        # Use the base class implementation
        #print(kwargs)
        curve = super(Relation, cls).from_columns(*columns, **kwargs)

        # Set x name and y name
        curve.x_name = curve.column_names[0]
        curve.y_name = curve.column_names[1]

        # Return the curve
        return curve

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Load the curve using the base class implementation
        relation = super(Relation, cls).from_file(path, **kwargs)

        # Set x name and y name
        relation.x_name = relation.column_info[0][0]
        relation.y_name = relation.column_info[1][0]

        # Return the relation
        return relation

    # -----------------------------------------------------------------

    @classmethod
    def from_data_file(cls, path, x_name="x", x_description="x values", x_unit=None, y_name="y",
                       y_description="y values", y_unit=None, x_column=0, y_column=1, skiprows=0, conversion_info=None):

        """
        This function ...
        :param path:
        :param x_name:
        :param x_description:
        :param x_unit:
        :param y_name:
        :param y_description:
        :param y_unit:
        :param x_column:
        :param y_column:
        :param skiprows:
        :param conversion_info:
        :return:
        """

        # Define columns
        columns = (x_column, y_column)

        kwargs = dict()

        kwargs["x_name"] = x_name
        kwargs["x_description"] = x_description
        kwargs["x_unit"] = x_unit

        kwargs["y_name"] = y_name
        kwargs["y_description"] = y_description
        kwargs["y_unit"] = y_unit

        # Create the curve
        curve = cls(**kwargs)
        x_values, y_values = np.loadtxt(path, unpack=True, usecols=columns, skiprows=skiprows)
        for x_value, y_value in zip(x_values, y_values): curve.add_point(x_value, y_value, conversion_info=conversion_info)

        # Return the curve
        return curve

    # -----------------------------------------------------------------

    @property
    def x_unit(self):
        return self[self.x_name].unit

    # -----------------------------------------------------------------

    @property
    def y_unit(self):
        return self[self.y_name].unit

    # -----------------------------------------------------------------

    @property
    def x_data(self):
        return self[self.x_name]

    # -----------------------------------------------------------------

    @property
    def x_array(self):
        return np.asarray(self.x_data)

    # -----------------------------------------------------------------

    @property
    def y_data(self):
        return self[self.y_name]

    # -----------------------------------------------------------------

    @property
    def y_array(self):
        return np.asarray(self.y_data)

    # -----------------------------------------------------------------

    @property
    def has_aux(self):
        if self.has_errors: return self.ncolumns > 4
        else: return self.ncolumns > 2

    # -----------------------------------------------------------------

    @property
    def aux_indices(self):
        from ..tools import sequences
        if self.has_aux: return sequences.indices_not_in(self.column_names, [self.x_name, self.y_name, errormin_name, errorplus_name])
        else: return []

    # -----------------------------------------------------------------

    @property
    def aux_names(self):
        return [self.column_names[index] for index in self.aux_indices]

    # -----------------------------------------------------------------

    @property
    def aux_units(self):
        units = self.column_units
        return [units[index] for index in self.aux_indices]

    # -----------------------------------------------------------------

    @property
    def naux(self):
        if self.has_errors: return self.ncolumns - 4
        else: return self.ncolumns - 2

    # -----------------------------------------------------------------

    def add_aux(self, name, values, unit=None, description=None, as_column=False):

        """
        This function ...
        :param name:
        :param values:
        :param unit:
        :param description:
        :param as_column:
        :return:
        """

        # Add column info
        if description is None: description = "auxilary axis"

        # Add
        self.add_col(name, values, dtype=float, unit=unit, description=description, as_column=as_column)

    # -----------------------------------------------------------------

    def add_point(self, x_value, y_value, conversion_info=None, sort=False, aux=None):

        """
        This function ...
        :param x_value:
        :param y_value:
        :param conversion_info:
        :param sort: DEFAULT IS FALSE HERE
        :param aux:
        :return:
        """

        # Set values
        values = [x_value, y_value]

        # Add aux values
        if self.has_aux:
            for aux_name in self.aux_names:
                if aux is None or aux_name not in aux: aux_value = None
                else: aux_value = aux[aux_name]
                values.append(aux_value)

        # No aux columns but aux is passed: error
        elif aux is not None: raise ValueError("There are no aux columns")

        # Add a row
        self.add_row(values, conversion_info=conversion_info)

        # Sort the table by the x values
        if sort: self.sort(self.x_name)

    # -----------------------------------------------------------------

    @property
    def has_errors(self):
        return errormin_name in self.colnames and errorplus_name in self.colnames

    # -----------------------------------------------------------------

    def get_x(self, unit=None, asarray=False, add_unit=True, conversion_info=None, density=False, brightness=False):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :param conversion_info:
        :param density:
        :param brightness:
        :return:
        """

        # Create and return
        if asarray: return arrays.plain_array(self[self.x_name], unit=unit, array_unit=self.column_unit(self.x_name),
                                      conversion_info=conversion_info, density=density, brightness=brightness)
        else: return arrays.array_as_list(self[self.x_name], unit=unit, add_unit=add_unit,
                                        array_unit=self.column_unit(self.x_name), conversion_info=conversion_info,
                                        density=density, brightness=brightness)

    # -----------------------------------------------------------------

    def get_x_value(self, index, unit=None, add_unit=True, conversion_info=None, density=False, brightness=False):

        """
        This function ...
        :param index:
        :param unit:
        :param add_unit:
        :param conversion_info:
        :param density:
        :param brightness:
        :return:
        """

        from ..units.parsing import parse_unit as u
        if unit is not None: unit = u(unit, density=density, brightness=brightness)
        return self.get_value(self.x_name, index, unit=unit, add_unit=add_unit, conversion_info=conversion_info)

    # -----------------------------------------------------------------

    def get_min_x_value(self, unit=None, add_unit=True, conversion_info=None, density=False, brightness=False,
                        ignore_zero=False, ignore_negatives=False):

        """
        This function ...
        :param unit:
        :param add_unit:
        :param conversion_info:
        :param density:
        :param brightness:
        :param ignore_zero:
        :param ignore_negatives:
        :return:
        """

        # Get unit
        from ..units.parsing import parse_unit as u
        if unit is not None: unit = u(unit, density=density, brightness=brightness)
        else: unit = self.x_unit

        # Get values as array
        values = self.get_x(unit=unit, asarray=True, conversion_info=conversion_info)

        # Mask?
        if ignore_zero: values = values[values != 0]
        if ignore_negatives: values = values[values >= 0]

        # Get minimum
        min_value = np.nanmin(values)

        # Return
        if add_unit: return min_value * unit
        else: return min_value

    # -----------------------------------------------------------------

    def get_max_x_value(self, unit=None, add_unit=True, conversion_info=None, density=False, brightness=False, ignore_zero=False):

        """
        This function ...
        :param unit:
        :param add_unit:
        :param conversion_info:
        :param density:
        :param brightness:
        :param ignore_zero:
        :return:
        """

        # Get unit
        from ..units.parsing import parse_unit as u
        if unit is not None: unit = u(unit, density=density, brightness=brightness)
        else: unit = self.x_unit

        # Get values as array
        values = self.get_x(unit=unit, asarray=True, conversion_info=conversion_info)

        # Mask?
        if ignore_zero: values = values[values != 0]

        # Get maximum
        min_value = np.nanmax(values)

        # Return
        if add_unit: return min_value * unit
        else: return min_value

    # -----------------------------------------------------------------

    def get_y(self, unit=None, asarray=False, add_unit=True, conversion_info=None, density=False, brightness=False):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :param conversion_info:
        :param density:
        :param brightness:
        :return:
        """

        # Create and return
        if asarray: return arrays.plain_array(self[self.y_name], unit=unit, array_unit=self.column_unit(self.y_name),
                                      conversion_info=conversion_info, density=density, brightness=brightness)
        else: return arrays.array_as_list(self[self.y_name], unit=unit, add_unit=add_unit,
                                        array_unit=self.column_unit(self.y_name), conversion_info=conversion_info,
                                        density=density, brightness=brightness)

    # -----------------------------------------------------------------

    def get_y_value(self, index, unit=None, add_unit=True, conversion_info=None, density=False, brightness=False):

        """
        This function ...
        :param index:
        :param unit:
        :param add_unit:
        :param conversion_info:
        :param density:
        :param brightness:
        :return:
        """

        from ..units.parsing import parse_unit as u
        if unit is not None: unit = u(unit, density=density, brightness=brightness)
        return self.get_value(self.y_name, index, unit=unit, add_unit=add_unit, conversion_info=conversion_info)

    # -----------------------------------------------------------------

    def get_min_y_value(self, unit=None, add_unit=True, conversion_info=None, density=False, brightness=False,
                        ignore_zero=False, ignore_negatives=False):

        """
        This function ...
        :param unit:
        :param add_unit:
        :param conversion_info:
        :param density:
        :param brightness:
        :param ignore_zero:
        :param ignore_negatives:
        :return:
        """

        # Get unit
        from ..units.parsing import parse_unit as u
        if unit is not None: unit = u(unit, density=density, brightness=brightness)
        else: unit = self.y_unit

        # Get values as array
        values = self.get_y(unit=unit, asarray=True, conversion_info=conversion_info)

        # Mask?
        if ignore_zero: values = values[values!=0]
        if ignore_negatives: values = values[values>=0]

        # Get minimum
        min_value = np.nanmin(values)

        # Return
        if add_unit: return min_value * unit
        else: return min_value

    # -----------------------------------------------------------------

    def get_max_y_value(self, unit=None, add_unit=True, conversion_info=None, density=False, brightness=False, ignore_zero=False):

        """
        This function ...
        :param unit:
        :param add_unit:
        :param conversion_info:
        :param density:
        :param brightness:
        :param ignore_zero:
        :return:
        """

        # Get unit
        from ..units.parsing import parse_unit as u
        if unit is not None: unit = u(unit, density=density, brightness=brightness)
        else: unit = self.y_unit

        # Get values as array
        values = self.get_y(unit=unit, asarray=True, conversion_info=conversion_info)

        # Mask?
        if ignore_zero: values = values[values != 0]

        # Get maximum
        min_value = np.nanmax(values)

        # Return
        if add_unit: return min_value * unit
        else: return min_value

    # -----------------------------------------------------------------

    @property
    def npoints(self):
        return len(self)

    # -----------------------------------------------------------------

    def get_indices_right(self, x_min, include=True):

        """
        This function ...
        :param x_min:
        :param include:
        :return:
        """

        # Get the values
        x_values = self.get_x()

        # Initialize
        indices = []

        # Loop over the x values
        for index, value in enumerate(x_values):

            # Checks
            if x_min is not None:

                if include:
                    if value < x_min: continue
                else:
                    if value <= x_min: continue

            # Add the index
            indices.append(index)

        # Return the indices
        return indices

    # -----------------------------------------------------------------

    def get_indices_left(self, x_max, include=True):

        """
        This function ...
        :param x_max:
        :param include:
        :return:
        """

        # Get the values
        x_values = self.get_x()

        # Initialize
        indices = []

        # Loop over the x values
        for index, value in enumerate(x_values):

            # Checks
            if include:
                if value > x_max: continue
            else:
                if value >= x_max: continue

            # Add the index
            indices.append(index)

        # Return the indices
        return indices

    # -----------------------------------------------------------------

    def get_indices(self, x_min=None, x_max=None, include_min=True, include_max=True):

        """
        This function ...
        :param x_min:
        :param x_max:
        :param include_min:
        :param include_max:
        :return:
        """

        # No limits given
        if x_min is None and x_max is None: return list(range(self.npoints))

        # Get the values
        x_values = self.get_x()

        # Initialize
        indices = []

        # Loop over the values
        for index, value in enumerate(x_values):

            # Checks
            if x_min is not None:
                if include_min:
                    if value < x_min: continue
                else:
                    if value <= x_min: continue
            if x_max is not None:
                if include_max:
                    if value > x_max: continue
                else:
                    if value >= x_max: continue

            # Add the index
            indices.append(index)

        # Return the indices
        return indices

    # -----------------------------------------------------------------

    def splice_right(self, x_min, include=True):

        """
        This function splices at the given x value and returns the right part
        :param x_min:
        :param include:
        :return:
        """

        # Get the indices
        indices = self.get_indices_right(x_min, include=include)

        # Return
        return self.splice_indices(indices)

    # -----------------------------------------------------------------

    def splice_left(self, x_max, include=True):

        """
        This function splics the curve and the given x value and returns the left part
        :param x_max:
        :param include:
        :return:
        """

        # Get the indices
        indices = self.get_indices_left(x_max, include=include)

        # Return
        return self.splice_indices(indices)

    # -----------------------------------------------------------------

    def splice(self, x_min=None, x_max=None, include_min=True, include_max=True):

        """
        This function ...
        :param x_min:
        :param x_max:
        :param include_min:
        :param include_max:
        :return:
        """

        # Get the indices
        indices = self.get_indices(x_min, x_max, include_min=include_min, include_max=include_max) # don't name arguments because of re-definition of function in WavelengthCurve class

        # Return
        return self.splice_indices(indices)

    # -----------------------------------------------------------------

    def splice_indices(self, indices):

        """
        This function ...
        :param indices:
        :return:
        """

        # Set the x and y unit
        x_unit = self.x_unit
        y_unit = self.y_unit

        # Get the values
        x_values = [self.get_value(self.x_name, index, unit=x_unit, add_unit=False) for index in indices]
        y_values = [self.get_value(self.y_name, index, unit=y_unit, add_unit=False) for index in indices]

        # Set the names and units
        names = (self.x_name, self.y_name)
        units = (x_unit, y_unit)

        # Create new curve
        return self.__class__.from_columns(x_values, y_values, names=names, units=units)

    # -----------------------------------------------------------------

    def get_x_splice(self, x_min=None, x_max=None, include_min=True, include_max=True, asarray=False, return_indices=False):

        """
        This function ...
        :param x_min:
        :param x_max:
        :param include_min:
        :param include_max:
        :param asarray:
        :param return_indices:
        :return:
        """

        # Get the indices
        indices = self.get_indices(x_min, x_max, include_min=include_min, include_max=include_max)  # don't name arguments because of re-definition of function in WavelengthCurve class

        # Get the values
        x_values = [self.get_value(self.x_name, index, unit=self.x_unit, add_unit=False) for index in indices]

        # Return
        if asarray: x_values =  np.array(x_values)
        if return_indices: return x_values, indices
        else: return x_values

    # -----------------------------------------------------------------

    def get_y_splice(self, x_min=None, x_max=None, include_min=True, include_max=True, asarray=False, return_indices=False):

        """
        This function ...
        :param x_min:
        :param x_max:
        :param include_min:
        :param include_max:
        :param asarray:
        :param return_indices:
        :return:
        """

        # Get the indices
        indices = self.get_indices(x_min, x_max, include_min=include_min, include_max=include_max)  # don't name arguments because of re-definition of function in WavelengthCurve class

        # Get the values
        y_values = [self.get_value(self.y_name, index, unit=self.y_unit, add_unit=False) for index in indices]

        # Return
        if asarray: return np.array(y_values)
        if return_indices: return y_values, indices
        else: return y_values

    # -----------------------------------------------------------------

    def flatten_above(self, value, flatten_value=0., include=True):

        """
        This function ...
        :param value:
        :param flatten_value:
        :param include:
        :return:
        """

        from ..units.parsing import parse_quantity

        # Check value with unit
        if self.x_unit is not None:
            if hasattr(value, "unit"): pass
            elif isinstance(value, basestring): value = parse_quantity(value)
            else: raise ValueError("Unit of the value is not defined")
        elif hasattr(value, "unit") or isinstance(value, basestring): raise ValueError("Unit of the value is defined, but column unit is not")

        # Get the indices of the values to be flattened
        #indices = self.get_indices(x_min=value, include_min=include)
        indices = self.get_indices(value, None, include_min=include) # works also with derived class implementation

        # Set flatten value with unit
        if self.y_unit is not None:
            if hasattr(flatten_value, "unit"): pass
            elif flatten_value == 0.: flatten_value = flatten_value * self.y_unit
            else: raise ValueError("Unit of the flatten value is not defined")
        elif hasattr(flatten_value, "unit") or isinstance(flatten_value, basestring): raise ValueError("Unit of the flatten value is defined, but column unit is not")
        #print(flatten_value)

        # Flatten values
        for index in indices: self.set_value(self.y_name, index, flatten_value)

    # -----------------------------------------------------------------

    def flatten_below(self, value, flatten_value=0., include=True):

        """
        This function ...
        :param value:
        :param flatten_value:
        :param include:
        :return:
        """

        from ..units.parsing import parse_quantity

        # Check value with unit
        if self.x_unit is not None:
            if hasattr(value, "unit"): pass
            elif isinstance(value, basestring): value = parse_quantity(value)
            else: raise ValueError("Unit of the value is not defined")
        elif hasattr(value, "unit") or isinstance(value, basestring): raise ValueError("Unit of the value is defined, but column unit is not")

        # Get the indices of the values to be flattened
        #indices = self.get_indices(x_max=value, include_max=include)
        indices = self.get_indices(None, value, include_max=include) # works also with derived class implementation

        # Set flatten value with unit
        if self.y_unit is not None:
            if hasattr(flatten_value, "unit"): pass
            elif isinstance(flatten_value, basestring): flatten_value = parse_quantity(flatten_value)
            elif flatten_value == 0.: flatten_value = flatten_value * self.y_unit
            else: raise ValueError("Unit of the flatten value is not defined")
        elif hasattr(flatten_value, "unit") or isinstance(flatten_value, basestring): raise ValueError("Unit of the flatten value is defined, but column unit is not")
        #print(flatten_value)

        # Flatten values
        for index in indices: self.set_value(self.y_name, index, flatten_value)

    # -----------------------------------------------------------------

    def flattened_above(self, value, flatten_value=0., include=True):

        """
        This function ...
        :param value:
        :param flatten_value:
        :param include:
        :return:
        """

        # Make copy
        new = self.copy()
        new.flatten_above(value, flatten_value=flatten_value, include=include)
        return new

    # -----------------------------------------------------------------

    def flattened_below(self, value, flatten_value=0., include=True):

        """
        This function ...
        :param value:
        :param flatten_value:
        :param include:
        :return:
        """

        # Make copy
        new = self.copy()
        new.flatten_below(value, flatten_value=flatten_value, include=include)
        return new

# -----------------------------------------------------------------
