#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.units.quantity Contains the PhotometricQuantity class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import math
import numpy as np

# Import astronomical modules
from astropy.units import Quantity
from astropy.coordinates import Angle

# Import the relevant PTS classes and modules
from ..tools import types
from .parsing import parse_quantity, parse_unit, parse_photometric_unit
from .stringify import represent_quantity
from .unit import PhotometricUnit, divide_units, multiply_units, get_conversion_factor

# -----------------------------------------------------------------

def get_zero_quantity(unit, density=False, brightness=False, density_strict=False, brightness_strict=False):

    """
    This function ...
    :param unit:
    :param density:
    :param brightness:
    :param density_strict:
    :param brightness_strict:
    :return:
    """

    return 0.0 * parse_unit(unit, density=density, brightness=brightness, density_strict=density_strict, brightness_strict=brightness_strict)

# -----------------------------------------------------------------

def get_zero_photometric_quantity(unit, density=False, brightness=False, density_strict=False, brightness_strict=False):

    """
    This function ...
    :param unit:
    :param density:
    :param brightness:
    :param density_strict:
    :param brightness_strict:
    :return:
    """

    return 0.0 * parse_photometric_unit(unit, density=density, brightness=brightness, density_strict=density_strict, brightness_strict=brightness_strict)

# -----------------------------------------------------------------

def zero_angle(): return Angle(0.0, "deg")

# -----------------------------------------------------------------

def right_angle(): return Angle(90., "deg")

# -----------------------------------------------------------------

def zero_length_quantity(): return 0.0 * parse_unit("m")

# -----------------------------------------------------------------

# NO: DOESN'T MAKE SENSE! DIFFERENT SCALES EXIST!!
#zero_temperature_quantity = 0.0 * parse_unit("K")

# -----------------------------------------------------------------

def zero_mass_quantity(): return 0.0 * parse_unit("kg")

# -----------------------------------------------------------------

def zero_mass_density_quantity(): return 0.0 * parse_unit("kg/m3")

# -----------------------------------------------------------------

def is_scalar(value):

    """
    This function ...
    :param value:
    :return:
    """

    if types.is_real_type(value): return True
    elif types.is_quantity(value): return value.unit == ""
    else: raise ValueError("Value type not recognized: must be real value or quantity")

# -----------------------------------------------------------------

def get_scalar(value):

    """
    This function ...
    :param value:
    :return:
    """

    if types.is_real_type(value): return value
    elif types.is_quantity(value): return value.value
    else: raise ValueError("Value type not recognized: must be real value or quantity")

# -----------------------------------------------------------------

class PhotometricQuantity(Quantity):
    
    """
    This class ...
    """
    
    def __new__(cls, value, unit):
        
        """
        The constructor ...
        :param value:
        :param unit:
        """

        from .unit import PhotometricUnit

        # Create unit
        if types.is_string_type(unit): unit = PhotometricUnit(unit)

        # Call constructor of the base class
        quantity = Quantity.__new__(cls, value, unit)

        # Set photometric unit, the base class (Quantity) converts the unit back to an ordinary Unit in the __new__ function above
        quantity._unit = unit

        # Return the new quantity
        return quantity

    # -----------------------------------------------------------------

    def __quantity_subclass__(self, unit):

        """
        Overridden by subclasses to change what kind of view is
        created based on the output unit of an operation.

        Parameters
        ----------
        unit : UnitBase
            The unit for which the appropriate class should be returned

        Returns
        -------
        tuple :
            - `Quantity` subclass
            - bool: True is subclasses of the given class are ok
        """

        #return Quantity, True

        from .parsing import parse_unit
        from .unit import PhotometricUnit

        unit = parse_unit(unit)
        if isinstance(unit, PhotometricUnit): return PhotometricQuantity, True
        else: return Quantity, True

    # -----------------------------------------------------------------

    @Quantity.unit.setter
    def unit(self, argument):

        """
        This function ...
        :param argument:
        :return:
        """

        from .unit import PhotometricUnit

        # Parse as photometric unit
        unit = PhotometricUnit(argument)
        self._unit = unit

    # -----------------------------------------------------------------

    @property
    def scale_factor(self):

        """
        This function ...
        :return:
        """

        return self.unit.scale_factor

    # -----------------------------------------------------------------

    def __add__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        # Call the implementation
        new_value = add_with_units(self.value, self.unit, other)

        # Create new quantity
        return PhotometricQuantity(new_value, self.unit)

    # -----------------------------------------------------------------

    def __sub__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        # Call the implementation
        new_value = subtract_with_units(self.value, self.unit, other)

        # Create new quantity
        return PhotometricQuantity(new_value, self.unit)

    # -----------------------------------------------------------------

    def __mul__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        # Call the implementation
        new_value, new_unit = multiply_with_units(self.value, self.unit, other)
        if new_unit is None: new_unit = parse_unit("") # dimensionless

        # Create new quantity
        return parse_quantity(new_value * new_unit)

    # -----------------------------------------------------------------

    def __div__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        # Call the implementation
        new_value, new_unit = divide_with_units(self.value, self.unit, other)
        if new_unit is None: new_unit = parse_unit("") # dimensionless

        # Create new quantity
        return parse_quantity(new_value * new_unit)

    # -----------------------------------------------------------------

    __truediv__ = __div__

    # -----------------------------------------------------------------

    @property
    def base_unit(self):
        return self.unit.base_unit

    # -----------------------------------------------------------------

    @property
    def wavelength_unit(self):
        return self.unit.wavelength_unit

    # -----------------------------------------------------------------

    @property
    def frequency_unit(self):
        return self.unit.frequency_unit

    # -----------------------------------------------------------------

    @property
    def distance_unit(self):
        return self.unit.distance_unit

    # -----------------------------------------------------------------

    @property
    def extent_unit(self):
        return self.unit.extent_unit

    # -----------------------------------------------------------------

    @property
    def solid_angle_unit(self):
        return self.unit.solid_angle_unit

    # -----------------------------------------------------------------

    def copy(self):

        """
        This function ...
        :return:
        """

        return PhotometricQuantity(self.value, self.unit.copy())

    # -----------------------------------------------------------------

    @property
    def physical_type(self):
        return self.unit.physical_type

    # -----------------------------------------------------------------

    @property
    def base_physical_type(self):
        return self.unit.base_physical_type

    # -----------------------------------------------------------------

    @property
    def base_symbol(self):
        return self.unit.base_symbol

    # -----------------------------------------------------------------

    @property
    def symbol(self):
        return self.unit.symbol

    # -----------------------------------------------------------------

    @property
    def density(self):
        return self.unit.density

    # -----------------------------------------------------------------

    @property
    def brightness(self):
        return self.unit.brightness

    # -----------------------------------------------------------------

    @property
    def spectral_density_type(self):
        return self.unit.spectral_density_type

    # -----------------------------------------------------------------

    @property
    def is_spectral_density(self):
        return self.unit.is_spectral_density

    # -----------------------------------------------------------------

    @property
    def is_bolometric(self):
        return self.unit.is_bolometric

    # -----------------------------------------------------------------

    @property
    def is_wavelength_density(self):
        return self.unit.is_wavelength_density

    # -----------------------------------------------------------------

    @property
    def is_frequency_density(self):
        return self.unit.is_frequency_density

    # -----------------------------------------------------------------

    @property
    def is_neutral_density(self):
        return self.unit.is_neutral_density

    # -----------------------------------------------------------------

    def to(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        # Get properties
        unit = args[0]
        density = kwargs.pop("density", False)
        density_strict = kwargs.pop("density_strict", False)
        brightness = kwargs.pop("brightness", False)
        brightness_strict = kwargs.pop("brightness_strict", False)
        wavelength = kwargs.pop("wavelength", None)
        frequency = kwargs.pop("frequency", None)
        distance = kwargs.pop("distance", None)
        solid_angle = kwargs.pop("solid_angle", None)
        fltr = kwargs.pop("fltr", kwargs.pop("filter", None))
        pixelscale = kwargs.pop("pixelscale", None)
        silent = kwargs.pop("silent", False)

        # Equivalencies
        equivalencies = kwargs.pop("equivalencies", None)

        # If equivalencies are specified, use the implementation of the base class
        if equivalencies is not None:

            # Re-evaluate whether it is a photometric quantity or not by passing a stringified version to the parser
            return parse_quantity(represent_quantity(Quantity.to(self, unit, equivalencies=equivalencies)))

        # Import
        from .unit import PhotometricUnit

        # Parse the new photometric unit
        unit = PhotometricUnit(unit, density=density, brightness=brightness, density_strict=density_strict, brightness_strict=brightness_strict)

        # Determine conversion factor
        factor = self.unit.conversion_factor(unit, wavelength=wavelength, frequency=frequency, distance=distance, solid_angle=solid_angle, fltr=fltr, pixelscale=pixelscale, silent=silent)

        # Return new quantity
        return PhotometricQuantity(self.value * factor, unit)

# -----------------------------------------------------------------

def quantity_to_angle(qty, unit="deg"):

    """
    This function ...
    :param qty:
    :param unit:
    :return:
    """

    return Angle(qty.to(unit).value, unit)

# -----------------------------------------------------------------

def same_physical_type(qty_a, qty_b):

    """
    This function ...
    :param qty_a:
    :param qty_b:
    :return:
    """

    return qty_a.unit.physical_type == qty_b.unit.physical_type

# -----------------------------------------------------------------

def add_with_units(value, unit, other, other_unit=None, conversion_info=None):

    """
    THis function ...
    :param value:
    :param unit:
    :param other:
    :param other_unit:
    :param conversion_info:
    :return:
    """

    # Does the first value have a unit?
    has_unit = unit is not None

    # Initialize conversion info
    if conversion_info is None: conversion_info = {}

    # Regular number, but unit defined?
    if types.is_real_or_integer(other) and other_unit is not None:
        other = other * other_unit
        other_unit = None

    # Regular number
    if types.is_real_or_integer(other):
        if has_unit: raise ValueError("First value has a unit")
        new_value = value + other

    # Quantity
    elif types.is_quantity(other):

        if other_unit is not None: raise ValueError("Cannot specify unit of second value when it is already a quantity")
        if not has_unit: raise ValueError("First value has no unit")
        new_value = value + other.to(unit, **conversion_info).value

    # Frame
    elif types.is_real_or_integer_array(other):

        # With unit
        if other_unit is not None:

            if not has_unit: raise ValueError("First value has no unit")
            #conversion_factor = other_unit.conversion_factor(unit, **conversion_info)
            # from_unit, to_unit, distance=None, wavelength=None, solid_angle=None, silent=False, parse=True, conversion_info=None
            conversion_factor = get_conversion_factor(other_unit, unit, conversion_info=conversion_info)
            new_value = value + conversion_factor * other

        # Without unit
        else:
            if has_unit: raise ValueError("First value has a unit")
            new_value = value + other

    # Else:
    else: raise ValueError("Invalid argument of type '" + str(type(other)) + "'")

    # Return the new value
    return new_value

# -----------------------------------------------------------------

def subtract_with_units(value, unit, other, other_unit=None, conversion_info=None):

    """
    This function ...
    :param value:
    :param unit:
    :param other:
    :param other_unit:
    :param conversion_info:
    :return:
    """

    # Does the first value have a unit?
    has_unit = unit is not None

    # Initialize conversion info
    if conversion_info is None: conversion_info = {}

    # Regular number, but unit defined?
    if types.is_real_or_integer(other) and other_unit is not None:
        other = other * other_unit # becomes quantity
        other_unit = None

    # Regular number
    if types.is_real_or_integer(other):
        if has_unit: raise ValueError("First value has a unit")
        new_value = value - other

    # Quantity
    elif types.is_quantity(other):

        if other_unit is not None: raise ValueError("Cannot specify unit of second value when it is already a quantity")
        if not has_unit: raise ValueError("First value has no unit")
        new_value = value - other.to(unit, **conversion_info).value

    # Frame
    elif types.is_real_or_integer_array(other):

        # With unit
        if other_unit is not None:

            if not has_unit: raise ValueError("First value has no unit")
            #self._data -= value.converted_to(self.unit).data
            #conversion_factor = other_unit.conversion_factor(unit, **conversion_info)
            conversion_factor = get_conversion_factor(other_unit, unit, conversion_info=conversion_info)
            new_value = value - conversion_factor * other

        # Without unit
        else:
            if has_unit: raise ValueError("First value has a unit")
            #self._data -= value.data
            new_value = value - other

    # Else
    else: raise ValueError("Invalid argument of type '" + str(type(other)) + "'")

    # Return the new value
    return new_value

# -----------------------------------------------------------------

def multiply_with_units(value, unit, other, other_unit=None):

    """
    This function ...
    :param value:
    :param unit:
    :param other:
    :param other_unit:
    :return:
    """

    has_unit = unit is not None

    # Regular number, but unit defined?
    if types.is_real_or_integer(other) and other_unit is not None:
        other = other * other_unit  # becomes quantity
        other_unit = None

    # Multiply by regular value
    if types.is_real_or_integer(other):

        # Determine value and unit
        new_value = value * float(other)
        new_unit = unit

    # Multiply by unit
    elif types.is_unit(other):

        other = parse_unit(other)
        if has_unit: new_unit = multiply_units(unit, other)
        else: new_unit = other

        # Photometric unit?
        if isinstance(new_unit, PhotometricUnit):

            if new_unit.has_scale:
                new_value = value * new_unit.scale_factor
                new_unit = new_unit.reduced_root  # without scale
            else:
                new_value = value
                new_unit = new_unit.reduced

        # Not photometric unit
        elif types.is_unit(new_unit):

            new_value = value * new_unit.scale
            new_unit._scale = 1

        # NOT A UNIT: UNITS MULTIPLIED AWAY
        elif types.is_real_or_integer(new_unit):

            new_value = value * new_unit
            new_unit = None

        # Invalid
        else: raise RuntimeError("Something went wrong")

    # Multiply by quantity
    elif types.is_quantity(other):

        if other_unit is not None: raise ValueError("Cannot specify unit of second value when it is already a quantity")

        other = parse_quantity(other)
        if has_unit: new_unit = multiply_units(unit, other.unit)
        else: new_unit = other.unit

        # Photometric unit?
        if isinstance(new_unit, PhotometricUnit):

            if new_unit.has_scale:
                new_value = value * other.value * new_unit.scale_factor
                new_unit = new_unit.reduced_root  # without scale
            else:
                new_value = value * other.value
                new_unit = new_unit.reduced

        # Not photometric unit
        elif types.is_unit(new_unit):

            new_value = value * other.value * new_unit.scale
            new_unit._scale = 1

        # NOT A UNIT: UNITS MULTPLIED AWAY (FOR EXAMPLE Hz * s)
        elif types.is_real_or_integer(new_unit):

            new_value = value * other.value * new_unit
            new_unit = None

        # Invalid
        else: raise RuntimeError("Something went wrong")

    # Multiply by array
    elif types.is_real_or_integer_array(other):

        # Both with unit
        if has_unit and other_unit is not None: new_unit = multiply_units(unit, other_unit)

        # Only first
        elif has_unit: new_unit = unit

        # Only other
        elif other_unit is not None: new_unit = other_unit

        # None
        else: new_unit = None

        # Photometric unit?
        if isinstance(new_unit, PhotometricUnit):

            if new_unit.has_scale:
                new_value = value * other * new_unit.scale_factor
                new_unit = new_unit.reduced_root  # without scale
            else:
                new_value = value * other
                new_unit = new_unit.reduced

        # Not photometric unit
        elif types.is_unit(new_unit):

            new_value = value * other.value * new_unit.scale
            new_unit._scale = 1

        # NOT A UNIT: UNITS MULTIPLIED AWAY (FOR EXAMPLE HZ * s)
        elif types.is_real_or_integer(new_unit):

            new_value = value * other.value * new_unit
            new_unit = None

        # Invalid
        else: raise RuntimeError("Something went wrong")

    # Invalid
    else: raise ValueError("Invalid argument of type '" + str(type(other)) + "'")

    # Return the new value and the new unit
    return new_value, new_unit

# -----------------------------------------------------------------

def get_dimension(data):

    """
    This function ...
    :param data:
    :return:
    """

    if types.is_datacube(data): return 3
    elif types.is_frame(data): return 2
    elif types.is_real_or_integer(data): return 1
    elif types.is_array_like(data): return data.ndim
    else: raise ValueError("Unknown type")

# -----------------------------------------------------------------

def float_division(a, b):

    """
    This function ...
    :param a:
    :param b:
    :return:
    """

    # Divide by zero, and not array like?
    if b == 0 and get_dimension(a) == 1: return math.copysign(float("inf"), a)

    # Create float of a
    if types.is_frame(a): pass
    elif types.is_array_like(a): a = a.astype(float)
    elif types.is_real_or_integer(a): a = float(a) # make regular float
    else: raise ValueError("Type of 'a' not recognized")

    # Create float of b
    if types.is_frame(b): pass
    elif types.is_array_like(b): b = b.astype(float)
    elif types.is_real_or_integer(b): b = float(b) # make regular float
    else: raise ValueError("Type of 'b' not recognized")

    # Return the result
    return a / b

# -----------------------------------------------------------------

def divide_with_units(value, unit, other, other_unit=None):

    """
    This function returns the result from dividing a 1D or higher-dimensional set of values, corresponding to a certain unit,
    with anything else, where a unit can also be expected
    :param value:
    :param unit:
    :param other:
    :param other_unit:
    :return:
    """

    has_unit = unit is not None

    # Regular number, but unit defined?
    if types.is_real_or_integer(other) and other_unit is not None:
        other = other * other_unit  # becomes quantity
        other_unit = None

    # Divide by regular value
    if types.is_real_or_integer(other):

        # Determine value and unit
        new_value = float_division(value, other)
        new_unit = unit

    # Divide by unit
    elif types.is_unit(other):

        other = parse_unit(other)
        if has_unit: new_unit = divide_units(unit, other)
        else: new_unit = 1./other

        # Photometric unit?
        if isinstance(new_unit, PhotometricUnit):

            if new_unit.has_scale:
                new_value = value * new_unit.scale_factor
                new_unit = new_unit.reduced_root  # without scale
            else:
                new_value = value
                new_unit = new_unit.reduced

        # Not photometric unit
        elif types.is_unit(new_unit):

            new_value = value * new_unit.scale
            new_unit._scale = 1

        # NOT A UNIT: UNITS DIVIDED AWAY?
        elif types.is_real_or_integer(new_unit):

            new_value = value * new_unit
            new_unit = None

        # Invalid new unit
        else: raise RuntimeError("Something went wrong")

    # Divide by quantity
    elif types.is_quantity(other):

        if other_unit is not None: raise ValueError("Cannot specify unit of second value when it is already a quantity")

        #print("UNIT", unit)
        #print("OTHER", other)

        other = parse_quantity(other)
        if has_unit: new_unit = divide_units(unit, other.unit)
        else: new_unit = other.unit**-1

        #print("NU", new_unit)

        # Photometric?
        if isinstance(new_unit, PhotometricUnit):

            if new_unit.has_scale:
                #if types.is_real_or_integer(value): new_value = float_division(value, other.value) * new_unit.scale_factor
                #elif types.is_real_or_integer_array(value): new_value = value / float(other.value) * new_unit.scale_factor
                #else: raise RuntimeError("Something went wrong")
                new_value = float_division(value, other.value) * new_unit.scale_factor
                new_unit = new_unit.reduced_root  # without scale
            else:
                #if types.is_real_or_integer(value): new_value = float_division(value, other.value)
                #elif types.is_real_or_integer_array(value): new_value = value / float(other.value)
                #else: raise RuntimeError("Something went wrong")
                new_value = float_division(value, other.value)
                new_unit = new_unit.reduced

        # Not photometric, but still unit
        elif types.is_unit(new_unit):

            #if types.is_real_or_integer(value): new_value = float_division(value, other.value) * new_unit.scale
            #elif types.is_real_or_integer_array(value): new_value = value / float(other.value) * new_unit.scale
            #else: raise RuntimeError("Something went wrong")
            new_value = float_division(value, other.value) * new_unit.scale
            new_unit._scale = 1

        # NOT A UNIT: UNITS DIVIDED AWAY?
        elif types.is_real_or_integer(new_unit):

            #if types.is_real_or_integer(value): new_value = float_division(value, other.value) * new_unit
            #elif types.is_real_or_integer_array(value): new_value = value / float(other.value) * new_unit
            #else: raise RuntimeError("Something went wrong")
            new_value = float_division(value, other.value) * new_unit
            new_unit = None

        # Invalid new unit
        else: raise RuntimeError("Something went wrong: type of unit: '" + type(new_unit).__name__ + "'")

    # Divide by array
    elif types.is_real_or_integer_array(other):

        # Both with unit
        if has_unit and other_unit is not None: new_unit = divide_units(unit, other_unit)

        # Only first
        elif has_unit: new_unit = unit

        # Only other
        elif other_unit is not None: new_unit = other_unit

        # None
        else: new_unit = None

        # Photometric?
        if isinstance(new_unit, PhotometricUnit):

            if new_unit.has_scale:
                new_value = value / other * new_unit.scale_factor
                new_unit = new_unit.reduced_root  # without scale
            else:
                new_value = value / other
                new_unit = new_unit.reduced

        # Not photometric
        elif types.is_unit(new_unit):

            new_value = value / other * new_unit.scale
            new_unit._scale = 1

        # NOT A UNIT: UNITS DIVIDED AWAY?
        elif types.is_real_or_integer(new_unit):

            new_value = value / other * new_unit
            new_unit = None

        # No unit
        elif new_unit is None: new_value = value / other

        # Invalid new unit
        else: raise RuntimeError("Something went wrong")

    # Invalid
    else: raise ValueError("Invalid argument of type '" + str(type(other)) + "'")

    # Return the new value and the new unit
    return new_value, new_unit

# -----------------------------------------------------------------

def get_value_and_unit(value):

    """
    This function splits a quantity or Frame into the actual data (or value) and the unit
    For a scalar, it returns the value, and None for the unit
    :param value:
    :return:
    """

    # Has unit property?
    if hasattr(value, "unit"):

        # Get unit
        unit = getattr(value, "unit")

        # Get value
        if hasattr(value, "value"): value = value.value # for quantity
        elif hasattr(value, "data"): value = value.data # for Frame
        elif hasattr(value, "values"): value = value.values # for Data3D
        else: pass  # assume using value is OK ...

    # No unit
    else: unit = None

    # Return
    return value, unit

# -----------------------------------------------------------------
