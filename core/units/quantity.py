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

# Import astronomical modules
from astropy.units import Quantity
from astropy.coordinates import Angle

# Import the relevant PTS classes and modules
from ..tools import types
from .parsing import parse_quantity, parse_unit, parse_photometric_unit
from .stringify import represent_quantity
from .unit import PhotometricUnit, divide_units, multiply_units

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

    def __div__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        # Divide by regular value
        if types.is_real_type(other) or types.is_integer_type(other):
            super(PhotometricQuantity, self).__div__(other)

        # Divide by unit
        elif types.is_unit(other):

            other = parse_unit(other)
            new_unit = divide_units(self.unit, other)

            # Photometric?
            if isinstance(new_unit, PhotometricUnit):
                if new_unit.has_scale:
                    new_value = self.value * new_unit.scale_factor
                    new_unit = new_unit.reduced_root # without scale
                else:
                    new_value = self.value
                    new_unit = new_unit.reduced

            # Not photometric
            else:
                new_value = self.value * new_unit.scale
                new_unit._scale = 1

            # Create new quantity
            return parse_quantity(new_value * new_unit)

        # Divide by quantity
        elif types.is_quantity(other):

            # Previous implementation
            #other = parse_quantity(other)
            # if isinstance(other, PhotometricQuantity):
            #     #print(self, self.physical_type)
            #     #print(other, other.physical_type)
            #     other_value = other.value
            #     other_unit = other.unit
            #     new_value = self.value / other_value
            #     new_unit = self.unit / other_unit
            #     return new_value * new_unit
            # # Let Astropy handle it ...
            # else: return super(PhotometricQuantity, self).__div__(other)

            other = parse_quantity(other)
            new_unit = divide_units(self.unit, other.unit)

            # Photometric?
            if isinstance(new_unit, PhotometricUnit):
                if new_unit.has_scale:
                    new_value = self.value / other.value * new_unit.scale_factor
                    new_unit = new_unit.reduced_root  # without scale
                else:
                    new_value = self.value / other.value
                    new_unit = new_unit.reduced

            # Not photometric
            else:
                new_value = self.value / other.value * new_unit.scale
                new_unit._scale = 1

            # Create new quantity
            return parse_quantity(new_value * new_unit)

        # Invalid
        else: raise ValueError("Invalid argument of type '" + str(type(other)) + "'")

    # -----------------------------------------------------------------

    def __mul__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        # if isinstance(other, Quantity):
        #     other = parse_quantity(other)
        #     if isinstance(other, PhotometricQuantity):
        #         other_value = other.value
        #         other_unit = other.unit

        # Multiply by regular value
        if types.is_real_type(other) or types.is_integer_type(other):
            super(PhotometricQuantity, self).__div__(other)

        # Multiply by unit
        elif types.is_unit(other):

            other = parse_unit(other)
            new_unit = multiply_units(self.unit, other)
            #print(new_unit)

            # Photometric?
            if isinstance(new_unit, PhotometricUnit):
                if new_unit.has_scale:
                    #print(new_unit.scale_factor)
                    #print(new_unit.reduced_root)
                    new_value = self.value * new_unit.scale_factor
                    new_unit = new_unit.reduced_root  # without scale
                else:
                    new_value = self.value
                    new_unit = new_unit.reduced

            # Not photometric
            else:
                new_value = self.value * new_unit.scale
                new_unit._scale = 1

            # Create new quantity
            return parse_quantity(new_value * new_unit)

        # Multiply by quantity
        elif types.is_quantity(other):

            other = parse_quantity(other)
            new_unit = multiply_units(self.unit, other.unit)
            #print(new_unit)

            # Photometric?
            if isinstance(new_unit, PhotometricUnit):
                if new_unit.has_scale:
                    new_value = self.value * other.value * new_unit.scale_factor
                    new_unit = new_unit.reduced_root  # without scale
                else:
                    new_value = self.value * other.value
                    new_unit = new_unit.reduced

            # Not photometric
            else:
                new_value = self.value * other.value * new_unit.scale
                new_unit._scale = 1

            # Create new quantity
            return parse_quantity(new_value * new_unit)

        # Invalid
        else: raise ValueError("Invalid argument of type '" + str(type(other)) + "'")

    # -----------------------------------------------------------------

    @property
    def base_unit(self):

        """
        This function ...
        :return:
        """

        return self.unit.base_unit

    # -----------------------------------------------------------------

    @property
    def wavelength_unit(self):

        """
        This function ...
        :return:
        """

        return self.unit.wavelength_unit

    # -----------------------------------------------------------------

    @property
    def frequency_unit(self):

        """
        This function ...
        :return:
        """

        return self.unit.frequency_unit

    # -----------------------------------------------------------------

    @property
    def distance_unit(self):

        """
        This function ...
        :return:
        """

        return self.unit.distance_unit

    # -----------------------------------------------------------------

    @property
    def extent_unit(self):

        """
        This function ...
        :return:
        """

        return self.unit.extent_unit

    # -----------------------------------------------------------------

    @property
    def solid_angle_unit(self):

        """
        This function ...
        :return:
        """

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

        """
        This function ...
        :return:
        """

        return self.unit.physical_type

    # -----------------------------------------------------------------

    @property
    def base_physical_type(self):

        """
        This function ...
        :return: 
        """

        return self.unit.base_physical_type

    # -----------------------------------------------------------------

    @property
    def base_symbol(self):

        """
        This function ...
        :return: 
        """

        return self.unit.base_symbol

    # -----------------------------------------------------------------

    @property
    def symbol(self):

        """
        This function ...
        :return: 
        """

        return self.unit.symbol

    # -----------------------------------------------------------------

    @property
    def density(self):

        """
        This function ...
        :return:
        """

        return self.unit.density

    # -----------------------------------------------------------------

    @property
    def brightness(self):

        """
        This function ...
        :return:
        """

        return self.unit.brightness

    # -----------------------------------------------------------------

    @property
    def spectral_density_type(self):

        """
        Thisn function ...
        :return:
        """

        return self.unit.spectral_density_type

    # -----------------------------------------------------------------

    @property
    def is_spectral_density(self):

        """
        This function ...
        :return:
        """

        return self.unit.is_spectral_density

    # -----------------------------------------------------------------

    @property
    def is_bolometric(self):

        """
        Thins function ...
        :return:
        """

        return self.unit.is_bolometric

    # -----------------------------------------------------------------

    @property
    def is_wavelength_density(self):

        """
        This function ...
        :return:
        """

        return self.unit.is_wavelength_density

    # -----------------------------------------------------------------

    @property
    def is_frequency_density(self):

        """
        This function ...
        :return:
        """

        return self.unit.is_frequency_density

    # -----------------------------------------------------------------

    @property
    def is_neutral_density(self):

        """
        Thisj function ...
        :return:
        """

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
        factor = self.unit.conversion_factor(unit, wavelength=wavelength, frequency=frequency, distance=distance, solid_angle=solid_angle, fltr=fltr, pixelscale=pixelscale)

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
