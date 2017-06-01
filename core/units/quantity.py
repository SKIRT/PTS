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
        wavelength = kwargs.pop("wavelength", None)
        frequency = kwargs.pop("frequency", None)
        distance = kwargs.pop("distance", None)
        solid_angle = kwargs.pop("solid_angle", None)
        fltr = kwargs.pop("fltr", None)
        pixelscale = kwargs.pop("pixelscale", None)

        # Equivalencies
        equivalencies = kwargs.pop("equivalencies", None)

        # If equivalencies are specified, use the implementation of the base class
        if equivalencies is not None:

            # Re-evaluate whether it is a photometric quantity or not by passing a stringified version to the parser
            return parse_quantity(represent_quantity(Quantity.to(self, unit, equivalencies=equivalencies)))

        # Import
        from .unit import PhotometricUnit

        #print(unit)

        # Parse the new photometric unit
        unit = PhotometricUnit(unit, density=density)

        #print(unit)

        # Determine conversion factor
        factor = self.unit.conversion_factor(unit, density=density, wavelength=wavelength, frequency=frequency, distance=distance, solid_angle=solid_angle, fltr=fltr, pixelscale=pixelscale)

        # Return new quantity
        return PhotometricQuantity(self.value * factor, unit)

# -----------------------------------------------------------------
