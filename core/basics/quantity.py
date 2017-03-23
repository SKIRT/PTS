#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.quantity Contains the PhotometricQuantity class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.units import Quantity
from astropy.coordinates import Angle

# -----------------------------------------------------------------

def parse_quantity(argument, density=False, physical_type=None):

    """
    This function ...
    :param argument:
    :param density:
    :param physical_type:
    :return:
    """

    from .unit import parse_unit

    if isinstance(argument, Quantity):

        number = argument.value
        unit = argument.unit
        unit = parse_unit(unit, density=density)

    else:

        # NEW IMPLEMENTATION
        units = ""
        number = 1.0
        while argument:
            try:
                number = float(argument)
                break
            except ValueError:
                units = argument[-1:] + units
                argument = argument[:-1]
        if units == "": raise ValueError("Unit is not specified")

        unit = parse_unit(units.strip(), density=density)

    # Check physical type
    if physical_type is not None:
        if unit.physical_type != physical_type: raise ValueError("The quantity was not parsed as a quantity of '" + physical_type + "'")

    # Create quantity
    return number * unit

# -----------------------------------------------------------------

def parse_angle(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    quantity = parse_quantity(argument, physical_type="angle")
    return Angle(quantity.value, quantity.unit)

# -----------------------------------------------------------------

def stringify_quantity(quantity):

    """
    This function ...
    :param quantity:
    :return:
    """

    from .unit import stringify_unit

    # Stringify the unit
    unit_type, unit_string = stringify_unit(quantity.unit)

    # Determine quantity parsing type
    if unit_type == "photometric_unit": parsing_type = "photometric_quantity"
    elif unit_type == "photometric_density_unit": parsing_type = "photometric_density_quantity"
    elif unit_type == "unit": parsing_type = "quantity"
    else: raise ValueError("Unknown unit type: " + unit_type)

    # Return parsing type and stringified quantity
    return parsing_type, repr(quantity.value) + " " + unit_string

# -----------------------------------------------------------------

def represent_quantity(quantity, scientific=False, decimal_places=2):

    """
    This function ...
    :param quantity:
    :param scientific:
    :param decimal_places:
    :return:
    """

    from .unit import represent_unit
    from ..tools.stringify import str_from_real
    return str_from_real(quantity.value, scientific=scientific, decimal_places=decimal_places) + " " + represent_unit(quantity.unit)

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
        if isinstance(unit, basestring): unit = PhotometricUnit(unit)

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

        from .unit import parse_unit as u
        from .unit import PhotometricUnit

        unit = u(unit)
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

        #return copy.deepcopy(self)
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
        equivalencies = kwargs.pop("equivalencies", None)

        # If equivalencies are specified, use the implementation of the base class
        if equivalencies is not None:

            # Re-evaluate whether it is a photometric quantity or not by passing a stringified version to the parser
            return parse_quantity(represent_quantity(Quantity.to(self, unit, equivalencies=equivalencies)))

        # Import
        from .unit import PhotometricUnit

        # Parse the new photometric unit
        unit = PhotometricUnit(unit, density=density)

        #print(self, self.unit, unit)

        print("")

        # Determine conversion factor
        factor = self.unit.conversion_factor(unit, density=density, wavelength=wavelength, frequency=frequency, distance=distance, solid_angle=solid_angle, fltr=fltr, pixelscale=pixelscale)

        from astropy.units import spectral
        if wavelength is not None:
            #print("")
            print("wavelength", wavelength)
            print("frequency", wavelength.to("Hz", equivalencies=spectral()))
            print("unit", self.unit)
            print("to unit", unit)
            print("conversion factor", factor)

        #print(wavelength, factor)

        # Return new quantity
        return PhotometricQuantity(self.value * factor, unit)

# -----------------------------------------------------------------

def is_mass_density(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    qty = parse_quantity(argument)
    return qty.unit.physical_type == "mass density"

# -----------------------------------------------------------------

def is_length(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    qty = parse_quantity(argument)
    return qty.unit.physical_type == "length"

# -----------------------------------------------------------------

def is_mass(argument):

    """
    This function ...
    :param argument:
    :return:
    """

    qty = parse_quantity(argument)
    return qty.unit.physical_type == "mass"

# -----------------------------------------------------------------

def is_photometric_density(argument):

    """
    This function ...
    :return:
    """

    qty = parse_quantity(argument, density=True)
    return qty.unit.is_spectral_density

# -----------------------------------------------------------------

def is_angle(argument):

    """
    This function ...
    :return:
    """

    qty = parse_quantity(argument)
    return qty.unit.physical_type == "angle"

# -----------------------------------------------------------------

def is_temperature(argument):

    """
    This function ...
    :return:
    """

    qty = parse_quantity(argument)
    return qty.unit.physical_type == "temperature"

# -----------------------------------------------------------------
