#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.projection.data Contains the DataProjections class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ...magic.core.frame import Frame
from ...magic.basics.vector import PixelShape
from ...core.tools import numbers
from ...core.basics.range import RealRange
from .general import default_scale_heights, is_faceon, is_edgeon, get_faceon_projection, get_edgeon_projection
from ...core.tools.utils import lazyproperty, ForbiddenOperation
from ...core.tools import sequences
from ...core.basics.log import log
from ...core.basics.table import SmartTable
from ...core.tools import types
from ...magic.basics.pixelscale import PhysicalPixelscale
from ...magic.basics.vector import Extent
from ...core.tools.stringify import tostr
from ...core.tools.progress import Bar, BAR_EMPTY_CHAR, BAR_FILLED_CHAR
from ...core.units.unit import PhotometricUnit
from ...core.units.parsing import parse_unit
from ...magic.core.frame import AllZeroError
from ...core.units.quantity import get_value_and_unit, add_with_units, subtract_with_units, multiply_with_units, divide_with_units

# -----------------------------------------------------------------

# Number of photon packages
default_npackages = 5e7

# Instruments/orientations
earth_name = "earth"
faceon_name = "faceon"
edgeon_name = "edgeon"

# -----------------------------------------------------------------

x_coordinate_name = "x coordinates"
y_coordinate_name = "y coordinate"
z_coordinate_name = "z coordinate"
weight_name = "weight"

# -----------------------------------------------------------------

class Data3D(object):

    """
    This class represents 3D data THAT IS STATIC! (not like a table, which is modifiable)
    """

    def __init__(self, name, x, y, z, values, weights=None, length_unit=None, unit=None, description=None, distance=None,
                 wavelength=None, solid_angle=None):

        """
        The constructor ...
        :param x:
        :param y:
        :param z:
        :param values:
        :param weights:
        :param length_unit:
        :param unit:
        :param description:
        :param distance:
        :param wavelength:
        :param solid_angle:
        """

        # Set the name and description for the data
        self.name = name
        self.description = description

        # Set coordinates
        self.x = x
        self.y = y
        self.z = z

        # Set values
        self.values = values

        # Set weights
        self.weights = weights

        # Set units
        self.length_unit = parse_unit(length_unit)
        self.unit = parse_unit(unit)

        # Set conversion info
        self.distance = distance
        self.wavelength = wavelength
        self.solid_angle = solid_angle

        # Check sizes?
        self.check()

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, length_unit=None, unit=None):

        """
        This function ...
        :param path:
        :param length_unit:
        :param unit:
        :return:
        """

        # Debugging
        log.debug("Loading the 3D data in table format from '" + path + "' ...")

        # Read the table
        table = SmartTable.from_file(path)

        # Find the name of the variable
        standard_column_names = [x_coordinate_name.capitalize(), y_coordinate_name.capitalize(), z_coordinate_name.capitalize(), weight_name.capitalize()]
        column_name = sequences.get_single_other(table.column_names, standard_column_names, none="error", method="error")
        name = column_name.lower()

        # Get the units
        if length_unit is None: length_unit = table.get_column_unit(x_coordinate_name.capitalize())
        if unit is None: unit = table.get_column_unit(column_name)

        # Get the data
        x = table.get_column_array(x_coordinate_name.capitalize(), unit=length_unit)
        y = table.get_column_array(y_coordinate_name.capitalize(), unit=length_unit)
        z = table.get_column_array(z_coordinate_name.capitalize(), unit=length_unit)
        values = table.get_column_array(name, unit=unit)
        weights = table.get_column_array(weight_name.capitalize()) if weight_name.capitalize() in table.colnames else None

        # Get the description
        description = table.meta["description"] if "description" in table.meta else None

        # Get the conversion info
        distance = table.meta["distance"] if "distance" in table.meta else None
        wavelength = table.meta["wavelength"] if "wavelength" in table.meta else None
        solid_angle = table.meta["solid_angle"] if "solid_angle" in table.meta else None

        # Create and return
        return cls(name, x, y, z, values, weights=weights, length_unit=length_unit, unit=unit, description=description,
                   distance=distance, wavelength=wavelength, solid_angle=solid_angle)

    # -----------------------------------------------------------------

    @property
    def has_description(self):
        return self.description is not None

    # -----------------------------------------------------------------

    @property
    def has_length_unit(self):
        return self.length_unit is not None

    # -----------------------------------------------------------------

    @property
    def has_unit(self):
        return self.unit is not None

    # -----------------------------------------------------------------

    @property
    def nx(self):
        return len(self.x)

    # -----------------------------------------------------------------

    @property
    def ny(self):
        return len(self.y)

    # -----------------------------------------------------------------

    @property
    def nz(self):
        return len(self.z)

    # -----------------------------------------------------------------

    @property
    def nvalues(self):
        return len(self.values)

    # -----------------------------------------------------------------

    @property
    def has_weights(self):
        return self.weights is not None

    # -----------------------------------------------------------------

    @property
    def nweights(self):
        return len(self.weights) if self.has_weights else 0

    # -----------------------------------------------------------------

    def __len__(self):

        """
        This function ...
        :return:
        """

        return self.nvalues

    # -----------------------------------------------------------------

    def check(self):

        """
        This function ...
        :return:
        """

        # Check types
        if not types.is_real_array(self.x): raise ValueError("x data must be an array of real values")
        if not types.is_real_array(self.y): raise ValueError("y data must be an array of real values")
        if not types.is_real_array(self.z): raise ValueError("z data must be an array of real values")
        if not types.is_real_array(self.values): raise ValueError("values must be an array of real values")
        if self.has_weights and not types.is_real_array(self.weights): raise ValueError("weights must be an array of real values")

        # Check
        sizes = [self.nx, self.ny, self.nz, self.nvalues]
        if not sequences.all_equal(sizes): raise ValueError("Sizes of data arrays are not equal")

        # Check weights
        if self.has_weights and self.nvalues != self.nweights: raise ValueError("Number of weights must be equal to number of values")

    # -----------------------------------------------------------------

    @lazyproperty
    def nans(self):
        return np.isnan(self.values)

    # -----------------------------------------------------------------

    @lazyproperty
    def infs(self):
        return np.isinf(self.values)

    # -----------------------------------------------------------------

    @lazyproperty
    def invalid(self):
        return self.nans + self.infs

    # -----------------------------------------------------------------

    @lazyproperty
    def masked_x(self):
        return np.ma.MaskedArray(self.x, mask=self.invalid)

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_x(self):
        return self.masked_x.compressed()

    # -----------------------------------------------------------------

    @lazyproperty
    def masked_y(self):
        return np.ma.MaskedArray(self.y, mask=self.invalid)

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_y(self):
        return self.masked_y.compressed()

    # -----------------------------------------------------------------

    @lazyproperty
    def masked_z(self):
        return np.ma.MaskedArray(self.z, mask=self.invalid)

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_z(self):
        return self.masked_z.compressed()

    # -----------------------------------------------------------------

    @lazyproperty
    def masked_values(self):
        return np.ma.MaskedArray(self.values, mask=self.invalid)

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_values(self):
        return self.masked_values.compressed()

    # -----------------------------------------------------------------

    @lazyproperty
    def masked_weights(self):
        if not self.has_weights: return None
        return np.ma.MaskedArray(self.weights, mask=self.invalid)

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_weights(self):
        if not self.has_weights: return None
        return self.masked_weights.compressed()

    # -----------------------------------------------------------------

    @property
    def has_distance(self):
        return self.distance is not None

    # -----------------------------------------------------------------

    @property
    def has_wavelength(self):
        return self.wavelength is not None

    # -----------------------------------------------------------------

    @property
    def has_solid_angle(self):
        return self.solid_angle is not None

    # -----------------------------------------------------------------

    @property
    def conversion_info(self):

        """
        Thisfunction ...
        :return:
        """

        info = dict()
        if self.has_distance: info["distance"] = self.distance
        if self.has_wavelength: info["wavelength"] = self.wavelength
        if self.has_solid_angle: info["solid_angle"] = self.solid_angle
        return info

    # -----------------------------------------------------------------

    def sum(self, add_unit=False):

        """
        This function ...
        :param add_unit:
        :return:
        """

        # Check
        unit = conversion_factor = None

        # Calculate
        result = np.nansum(self.values)

        # Convert?
        if conversion_factor is not None: result *= conversion_factor

        # Set unit
        if unit is None: unit = self.unit

        # Return
        if add_unit and self.has_unit: return result * unit
        else: return result

    # -----------------------------------------------------------------

    def __add__(self, other):

        """
        This fucntion ...
        :param other:
        :return:
        """

        # Get data and unit of other
        other_value, other_unit = get_value_and_unit(other)

        # Get new data
        new_values = add_with_units(self.values, self.unit, other_value, other_unit=other_unit, conversion_info=self.conversion_info)

        # Copy weights
        weights = np.copy(self.weights) if self.has_weights else None

        # Return the resulting data
        return self.__class__(self.name, self.x, self.y, self.z, new_values, weights=weights, length_unit=self.length_unit,
                              unit=self.unit, description=self.description, distance=self.distance,
                              wavelength=self.wavelength, solid_angle=self.solid_angle)

    # -----------------------------------------------------------------

    def __iadd__(self, other):
        # STATIC!
        raise ForbiddenOperation(self.__class__, "addition")

    # -----------------------------------------------------------------

    def __sub__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        # Get data and unit of other
        other_value, other_unit = get_value_and_unit(other)

        # Get new data
        new_values = subtract_with_units(self.values, self.unit, other_value, other_unit=other_unit, conversion_info=self.conversion_info)

        # Copy weights
        weights = np.copy(self.weights) if self.has_weights else None

        # Return the resulting data
        return self.__class__(self.name, self.x, self.y, self.z, new_values, weights=weights, length_unit=self.length_unit,
                              unit=self.unit, description=self.description, distance=self.distance,
                              wavelength=self.wavelength, solid_angle=self.solid_angle)

    # -----------------------------------------------------------------

    def __isub__(self, other):
        # STATIC!
        raise ForbiddenOperation(self.__class__, "subtraction")

    # -----------------------------------------------------------------

    def __mul__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        # Get data and unit of other
        other_value, other_unit = get_value_and_unit(other)

        # Get new data
        new_values, new_unit = multiply_with_units(self.values, self.unit, other_value, other_unit=other_unit)

        # Copy weights
        weights = np.copy(self.weights) if self.has_weights else None

        # Return the resulting data
        return self.__class__(self.name, self.x, self.y, self.z, new_values, weights=weights, length_unit=self.length_unit,
                              unit=new_unit, description=self.description, distance=self.distance,
                              wavelength=self.wavelength, solid_angle=self.solid_angle)

    # -----------------------------------------------------------------

    def __imul__(self, other):
        # STATIC!
        raise ForbiddenOperation(self.__class__, "multiplication")

    # -----------------------------------------------------------------

    def __div__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        # Get data and unit of other
        other_value, other_unit = get_value_and_unit(other)

        # Get new data
        new_values, new_unit = divide_with_units(self.values, self.unit, other_value, other_unit=other_unit)

        # Copy weights
        weights = np.copy(self.weights) if self.has_weights else None

        # Return the resulting data
        return self.__class__(self.name, self.x, self.y, self.z, new_values, weights=weights, length_unit=self.length_unit,
                              unit=new_unit, description=self.description, distance=self.distance,
                              wavelength=self.wavelength, solid_angle=self.solid_angle)

    # -----------------------------------------------------------------

    def __idiv__(self, other):
        # STATIC!
        raise ForbiddenOperation(self.__class__, "division")

    # -----------------------------------------------------------------

    __truediv__ = __div__

    # -----------------------------------------------------------------

    __itruediv__ = __idiv__

    # -----------------------------------------------------------------

    def copy(self, copy_coordinates=False):

        """
        This function ...
        :param copy_coordinates:
        :return:
        """

        # Make copies of the data
        values = np.copy(self.values)
        weights = np.copy(self.weights) if self.has_weights else None
        x = np.copy(self.x) if copy_coordinates else self.x
        y = np.copy(self.y) if copy_coordinates else self.y
        z = np.copy(self.z) if copy_coordinates else self.z

        # Return the new copy
        return self.__class__(self.name, x, y, z, values, weights=weights, length_unit=self.length_unit, unit=self.unit,
                              description=self.description, distance=self.distance, wavelength=self.wavelength,
                              solid_angle=self.solid_angle)

    # -----------------------------------------------------------------

    def normalized(self, to=1., return_sum=False, silent=False):

        """
        This function ...
        :param to:
        :param return_sum:
        :param silent:
        :return:
        """

        # Get the sum
        sum = self.sum()

        # Check whether the sum is nonnegative
        if sum < 0: raise RuntimeError("The sum of the data is negative")

        # Check if the sum is not zero
        if sum == 0: raise AllZeroError("The data cannot be normalized")

        # Calculate the conversion factor
        if hasattr(to, "unit"):  # quantity
            factor = to.value / sum
            unit = to.unit
        else:
            factor = to / sum
            unit = None

        # Debugging
        if not silent: log.debug("Multiplying the data with a factor of " + tostr(factor) + " to normalize to " + tostr(to) + " ...")

        # Create the normalized data and return
        new = self.converted_by_factor(factor, unit)
        if return_sum: return new, sum
        else: return new

    # -----------------------------------------------------------------

    def converted_to(self, to_unit, distance=None, density=False, brightness=False, density_strict=False,
                     brightness_strict=False, wavelength=None, solid_angle=None, silent=False):

        """
        This function ...
        :param to_unit:
        :param distance:
        :param density:
        :param brightness:
        :param density_strict:
        :param brightness_strict:
        :param wavelength:
        :param solid_angle:
        :param silent:
        :return:
        """

        # Check the unit of the data is defined
        if not self.has_unit: raise ValueError("The unit of the data is not defined")

        # Parse "to unit": VERY IMPORTANT, BECAUSE DOING SELF.UNIT = TO_UNIT WILL OTHERWISE REPARSE AND WILL BE OFTEN INCORRECT!! (NO DENSITY OR BRIGHTNESS INFO)
        to_unit = parse_unit(to_unit, density=density, brightness=brightness, brightness_strict=brightness_strict, density_strict=density_strict)

        # Already in the correct unit
        if to_unit == self.unit:
            if not silent: log.debug("Data is already in the desired unit")
            return self.copy()

        # Debugging
        if not silent: log.debug("Converting the data from unit " + tostr(self.unit, add_physical_type=True) + " to unit " + tostr(to_unit, add_physical_type=True) + " ...")

        # Get the conversion factor
        factor = self._get_conversion_factor(to_unit, distance=distance, wavelength=wavelength, solid_angle=solid_angle, silent=silent)

        # Debugging
        if not silent: log.debug("Conversion factor: " + str(factor))

        # Return converted data
        return self.converted_by_factor(factor, to_unit)

    # -----------------------------------------------------------------

    def converted_by_factor(self, factor, new_unit):

        """
        This function ...
        :param factor:
        :param new_unit:
        :return:
        """

        # Create multiplicated data
        new = self * factor

        # Set the new unit
        new.unit = new_unit

        # Return the new data instance
        return new

    # -----------------------------------------------------------------

    @property
    def is_photometric(self):
        return self.has_unit and isinstance(self.unit, PhotometricUnit)

    # -----------------------------------------------------------------

    def _get_conversion_factor(self, to_unit, distance=None, wavelength=None, solid_angle=None, silent=False):

        """
        This function ...
        :param to_unit:
        :param distance:
        :param wavelength:
        :param solid_angle:
        :param silent:
        :return:
        """

        # The data has a photometric unit
        if self.is_photometric:

            # Check that the target unit is also photometric
            if not isinstance(to_unit, PhotometricUnit): raise ValueError("Target unit is not photometric, while the data is")

            # Set the conversion info
            #conversion_info = dict()
            conversion_info = self.conversion_info
            if distance is not None: conversion_info["distance"] = distance
            if wavelength is not None: conversion_info["wavelength"] = wavelength
            if solid_angle is not None: conversion_info["solid_angle"] = solid_angle

            # Calculate the conversion factor
            factor = self.unit.conversion_factor(to_unit, silent=silent, **conversion_info)

        # The data does not have a photometric unit
        else:

            # Check whether target unit is also not photometric
            if isinstance(to_unit, PhotometricUnit): raise ValueError("Target unit is photometric, while the data is not")

            # Calculate the conversion factor
            factor = self.unit.to(to_unit, silent=True)

        # Return
        return factor

    # -----------------------------------------------------------------

    def saveto(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Debugging
        log.debug("Saving the 3D data in table format to '" + path + "' ...")

        # Set the columns
        columns = [self.x, self.y, self.z, self.values]
        names = [x_coordinate_name, y_coordinate_name, z_coordinate_name, self.name]
        if self.has_weights:
            columns.append(self.weights)
            columns.append(weight_name)

        # Capitalize column names
        names = [name.capitalize() for name in names]

        # Set units
        units = [self.length_unit, self.length_unit, self.length_unit, self.unit]
        if self.has_weights: units.append(None)

        # Create table
        table = SmartTable.from_columns(*columns, names=names, units=units, as_columns=True)

        # Set the description
        if self.has_description: table.meta["description"] = self.description

        # Set the conversion info
        if self.has_distance: table.meta["distance"] = self.distance
        if self.has_wavelength: table.meta["wavelength"] = self.wavelength
        if self.has_solid_angle: table.meta["solid_angle"] = self.solid_angle

        # Save the table
        table.saveto(path)

# -----------------------------------------------------------------

def project_3d(name, x_coordinates, y_coordinates, z_coordinates, values, projection, length_unit, unit=None,
               return_stddev=False, return_ncells=False, weights=None, description=None):

    """
    This function projects 3D data into a 2D map
    :param name: name of the kind of data
    :param x_coordinates:
    :param y_coordinates:
    :param z_coordinates:
    :param values:
    :param projection:
    :param length_unit:
    :param unit:
    :param return_stddev:
    :param return_ncells:
    :param weights:
    :param description:
    :return:
    """

    faceon = is_faceon(projection)
    edgeon = is_edgeon(projection)
    if faceon:
        faceon_projection = get_faceon_projection(projection, strict=True)
        edgeon_projection = None
    elif edgeon:
        faceon_projection = None
        edgeon_projection = get_edgeon_projection(projection, strict=True)
    else: raise ValueError("Projection must be face-on or edge-on")

    # Create data
    data = Data3D(x_coordinates, y_coordinates, z_coordinates, values, weights=weights, length_unit=length_unit, unit=unit)

    # Create the projections object
    projections = DataProjections(name, data, description=description, faceon=faceon, edgeon=edgeon, projection_faceon=faceon_projection, projection_edgeon=edgeon_projection)

    # Get the projected map
    if faceon: frame = projections.faceon
    elif edgeon: frame = projections.edgeon
    else: raise RuntimeError("Something went wrong")

    # Get the stddev
    if return_stddev:
        if faceon: stddev_frame = projections.faceon_stddev
        elif edgeon: stddev_frame = projection.edgeon_stddev
        else: raise RuntimeError("Something went wrong")
    else: stddev_frame = None

    # Get the ncells
    if return_ncells:
        if faceon: ncells_frame = projections.faceon_ncells
        elif edgeon: ncells_frame = projections.edgeon_ncells
        else: raise RuntimeError("Something went wrong")
    else: ncells_frame = None

    # Return
    if return_stddev:
        if return_ncells: return frame, stddev_frame, ncells_frame
        else: return frame, stddev_frame
    elif return_ncells: return frame, ncells_frame
    else: return frame

# -----------------------------------------------------------------

class DataProjections(object):

    """
    This class ...
    """

    def __init__(self, name, data, projection_faceon=None, projection_edgeon=None,
                 faceon=True, edgeon=True, description=None, distance=None, wcs=None, center=None,
                 radial_factor=1, scale_heights=default_scale_heights):

        """
        The constructor ...
        :param name:
        :param data:
        :param projection_faceon:
        :param projection_edgeon:
        :param path:
        :param faceon:
        :param edgeon:
        :param description:
        :param distance:
        :param wcs:
        :param center:
        :param radial_factor:
        :param scale_heights:
        """

        # Set
        self.name = name
        self.data = data

        # Set other
        self.description = description

        # Other
        self.distance = distance

        # Set the projections
        self.projection_faceon = projection_faceon
        self.projection_edgeon = projection_edgeon

        # Create projections?
        if faceon and not self.has_projection_faceon: self.create_projection_faceon()
        if edgeon and not self.has_projection_edgeon: self.create_projection_edgeon()

        # Create?
        if faceon: self.project_faceon()
        if edgeon: self.project_edgeon()

        # The maps
        self.faceon = None
        self.edgeon = None

        # The stddev maps
        self.faceon_stddev = None
        self.edgeon_stddev = None

        # The ncells maps
        self.faceon_ncells = None
        self.edgeon_ncells = None

    # -----------------------------------------------------------------

    @property
    def unit(self):
        return self.data.unit

    # -----------------------------------------------------------------

    @property
    def length_unit(self):
        return self.data.length_unit

    # -----------------------------------------------------------------

    @property
    def has_projection_faceon(self):
        return self.projection_faceon is not None

    # -----------------------------------------------------------------

    @property
    def has_projection_edgeon(self):
        return self.projection_edgeon is not None

    # -----------------------------------------------------------------

    def create_projection_faceon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the face-on projection ...")

    # -----------------------------------------------------------------

    def create_projection_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the edge-on projection ...")

    # -----------------------------------------------------------------

    @property
    def faceon_nx(self):
        return self.projection_faceon.pixels_x

    # -----------------------------------------------------------------

    @property
    def faceon_ny(self):
        return self.projection_faceon.pixels_y

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_npixels(self):
        return self.faceon_nx * self.faceon_ny

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_shape(self):
        return PixelShape.from_xy(self.faceon_nx, self.faceon_ny)

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_x_min(self):
        return self.projection_faceon.center_x - 0.5 * self.projection_faceon.field_x

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_x_max(self):
        return self.projection_faceon.center_x + 0.5 * self.projection_faceon.field_y

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_y_min(self):
        return self.projection_faceon.center_y - 0.5 * self.projection_faceon.field_y

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_y_max(self):
        return self.projection_faceon.center_y + 0.5 * self.projection_faceon.field_y

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_x_min_value(self):
        return self.faceon_x_min.to(self.length_unit).value

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_x_max_value(self):
        return self.faceon_x_max.to(self.length_unit).value

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_y_min_value(self):
        return self.faceon_y_min.to(self.length_unit).value

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_y_max_value(self):
        return self.faceon_y_max.to(self.length_unit).value

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_pixelscale(self):
        return self.projection_faceon.physical_pixelscale

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_pixelscale_value(self):
        return Extent(self.faceon_pixelscale.x.to(self.length_unit).value, self.faceon_pixelscale.y.to(self.length_unit).value)

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_half_pixelscale(self):
        return self.faceon_pixelscale * 0.5

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_half_pixelscale_value(self):
        return Extent(self.faceon_half_pixelscale.x.to(self.length_unit).value, self.faceon_half_pixelscale.y.to(self.length_unit).value)

    # -----------------------------------------------------------------

    def get_faceon_coordinate_ranges(self, x, y):

        """
        This function ...
        :param x:
        :param y:
        :return:
        """

        # Determine range of x and y
        x_range = RealRange(x - self.faceon_half_pixelscale_value.x, x + self.faceon_half_pixelscale_value.x)
        y_range = RealRange(y - self.faceon_half_pixelscale_value.y, y + self.faceon_half_pixelscale_value.y)

        # Return
        return x_range, y_range

    # -----------------------------------------------------------------

    @property
    def faceon_distance(self):
        return self.projection_faceon.distance

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_metadata(self):

        """
        Thisf unction ...
        :return:
        """

        meta = OrderedDict()
        meta["x_min"] = repr(x_min.value) + " " + tostr(self.length_unit)
        meta["x_max"] = repr(x_max.value) + " " + tostr(self.length_unit)
        meta["y_min"] = repr(y_min.value) + " " + tostr(self.length_unit)
        meta["y_max"] = repr(y_max.value) + " " + tostr(self.length_unit)
        return meta

    # -----------------------------------------------------------------

    def get_faceon_coordinate_mask(self, x, y):

        """
        This function ...
        :param x:
        :param y:
        :return:
        """

        # Get range
        x_range, y_range = self.get_faceon_coordinate_ranges(x, y)

        # Create masks
        x_mask = (x_range.min < valid_x_coordinates) * (valid_x_coordinates <= x_range.max)
        y_mask = (y_range.min < valid_y_coordinates) * (valid_y_coordinates <= y_range.max)

        # x_mask = self.get_coordinate_mask_x_for_map(x_range)
        # y_mask = self.get_coordinate_mask_y_for_map(y_range)
        # return x_mask * y_mask
        return x_mask * y_mask

    # -----------------------------------------------------------------

    def get_faceon_coordinate_indices(self, x, y):

        """
        This function ...
        :param x:
        :param y:
        :return:
        """

        # Get mask
        mask = self.get_faceon_coordinate_mask(x, y)

        # Get and return
        return np.where(mask)[0]

    # -----------------------------------------------------------------

    def project_faceon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the face-on map ...")

        # Initialize maps
        frame = Frame.initialize_nans(self.faceon_shape, unit=self.unit)
        stddev_frame = Frame.initialize_nans(self.faceon_shape, unit=self.unit)
        ncells_frame = Frame.initialize_nans(self.faceon_shape)

        # Set the pixelscale and the coordinate info
        frame.pixelscale = self.faceon_pixelscale
        frame.distance = self.faceon_distance
        frame.metadata.update(self.faceon_metadata)

        # Show progress bar
        with Bar(label='', width=32, hide=None, empty_char=BAR_EMPTY_CHAR,
                 filled_char=BAR_FILLED_CHAR, expected_size=self.faceon_npixels, every=1, add_datetime=True) as bar:

            # Loop over the pixels of the map
            x = self.faceon_x_min_value
            y = self.faceon_y_min_value
            index = 0
            for i, j in sequences.multirange(self.faceon_nx, self.faceon_ny):

                # Debugging
                # if index % 100 == 0: log.debug("Calculating heating fraction in the pixel " + str(index) + " of " + str(self.map_npixels) + " (" + tostr(float(index) / self.map_npixels * 100, decimal_places=1, round=True) + "%) ...")

                # Show progress
                progress = int(float(index+1) / float(self.faceon_npixels))
                bar.show(progress)

                # Get the indices
                indices = self.get_faceon_coordinate_indices(x, y)
                nindices = indices.shape[0]

                # Set number of cells
                ncells_frame[j, i] = nindices

                # If any cells
                if nindices > 0:

                    # Calculate the heating fraction
                    # fractions = self.valid_heating_fractions[indices]
                    # weights = self.valid_cell_weights[indices]
                    vals = valid_values[indices]
                    wghts = valid_weights[indices] if valid_weights is not None else None

                    # Calculate the mean heating fraction
                    fraction = numbers.weighed_arithmetic_mean_numpy(vals, weights=wghts)
                    fraction_stddev = numbers.weighed_standard_deviation_numpy(vals, weights=wghts, mean=fraction)

                    # Set fraction
                    frame[j, i] = fraction
                    stddev_frame[j, i] = fraction_stddev

                # Increment the x and y coordinate with one pixelsize
                x += self.faceon_pixelscale.x
                y += self.faceon_pixelscale.y
                index += 1

    # -----------------------------------------------------------------

    @property
    def edgeon_nx(self):
        return self.projection_edgeon.pixels_x

    # -----------------------------------------------------------------

    @property
    def edgeon_ny(self):
        return self.projection_edgeon.pixels_y

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_shape(self):
        return PixelShape.from_xy(self.edgeon_nx, self.edgeon_ny)

    # -----------------------------------------------------------------

    def project_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the edge-on map ...")

# -----------------------------------------------------------------
