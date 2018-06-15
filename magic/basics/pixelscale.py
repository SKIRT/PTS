#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.basics.pixelscale Contains the Pixelscale class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
from astropy.units import Unit, dimensionless_angles
from astropy.coordinates import Angle
from astropy.units.core import UnitTypeError, UnitConversionError

# Import standard modules
from .stretch import AngleExtent, PhysicalExtent

# -----------------------------------------------------------------

def angular_or_physical_pixelscale(*args, **kwargs):

    """
    This function ...
    :param args:
    :param kwargs:
    :return:
    """

    if len(args) == 1:
        if isinstance(args[0], Pixelscale): return args[0]
        elif isinstance(args[0], PhysicalPixelscale): return args[0]

    from ...core.units.parsing import parse_quantity
    from ...core.tools import types

    # Too many arguments
    if len(args) > 2: raise ValueError("Too many arguments")

    # Parse
    x = args[0] if len(args) > 0 else kwargs.get("x")
    y = args[1] if len(args) > 1 else kwargs.get("y", None)
    #print(x, y)
    x = parse_quantity(x, allow_scalar=True)
    y = parse_quantity(y, allow_scalar=True) if y is not None else None

    # Return pixelscale
    if types.is_length_quantity(x): return PhysicalPixelscale(x, y=y, **kwargs)
    elif types.is_angle(x): return Pixelscale(x, y=y, **kwargs)
    else:

        #print(x, y, kwargs)

        # Angular pixelscale
        try: pixelscale = Pixelscale(x, y=y, **kwargs)

        # Not angular pixelscale
        except (UnitTypeError, UnitConversionError) as e:

            # Physical pixelscale
            try: pixelscale = PhysicalPixelscale(x, y=y, **kwargs)

            # Invalid
            except (UnitTypeError, UnitConversionError) as e: raise ValueError("The passed value does not represent an angular or physical pixelscale")

        # Return
        return pixelscale

# -----------------------------------------------------------------

class Pixelscale(AngleExtent):

    """
    This class ...
    """

    def __init__(self, x, y=None, unit=None):

        """
        The constructor ...
        """

        # Remove '/pix' if present
        if y is None: y = x

        # Convert to just arcsec or just
        x = only_angle(x, unit=unit)
        y = only_angle(y, unit=unit)

        # Call the constructor of the base class
        super(Pixelscale, self).__init__(x, y)

    # -----------------------------------------------------------------

    @classmethod
    def from_physical(cls, pixelscale, distance, unit="arcsec"):

        """
        This function ...
        :param pixelscale:
        :param distance:
        :param unit:
        :return:
        """

        return cls(pixelscale.x_angle(distance, unit=unit), pixelscale.y_angle(distance, unit=unit))

    # -----------------------------------------------------------------

    @property
    def abs_x(self):

        """
        This function ...
        :return:
        """

        return abs(self.x)

    # -----------------------------------------------------------------

    @property
    def abs_y(self):

        """
        This function ...
        :return:
        """

        return abs(self.y)

    # -----------------------------------------------------------------

    @property
    def average(self):

        """
        This function ...
        :return:
        """

        x_pixelscale = abs(self.x.to("arcsec"))
        y_pixelscale = abs(self.y.to("arcsec"))

        # Check if the x pixelscale and y pixelscale are close
        if not np.isclose(x_pixelscale.value, y_pixelscale.value, rtol=0.001):

            from ...core.tools.stringify import tostr

            # Warn about the difference in x and y pixelscale
            from ...core.basics.log import log
            log.warning("Averaging the pixelscale over the x and y direction may not be a good approximation:")
            log.warning("  * x pixelscale (absolute value) = " + tostr(x_pixelscale))
            log.warning("  * y pixelscale (absolute value) = " + tostr(y_pixelscale))

        # Return a single value for the pixelscale in arcseconds
        return 0.5 * (x_pixelscale + y_pixelscale)

        # Other way from wcs:
        #np.mean(np.abs(np.diagonal(img_wcs.pixel_scale_matrix)))

    # -----------------------------------------------------------------

    @property
    def solid_angle(self):

        """
        This function ...
        :return:
        """

        solid_angle = (self.abs_x * self.abs_y).to("sr")
        return solid_angle

    # -----------------------------------------------------------------

    def x_extent(self, distance, unit="kpc"):

        """
        This function ...
        :param distance:
        :param unit:
        :return: 
        """

        length_x = (self.abs_x * distance).to(unit, equivalencies=dimensionless_angles())
        return length_x

    # -----------------------------------------------------------------

    def y_extent(self, distance, unit="kpc"):

        """
        This function ...
        :param distance:
        :param unit:
        :return: 
        """

        length_y = (self.abs_y * distance).to(unit, equivalencies=dimensionless_angles())
        return length_y

    # -----------------------------------------------------------------

    def extent(self, distance):

        """
        This function ...
        :param distance: 
        :return: 
        """

        from .stretch import PhysicalStretch
        return PhysicalStretch(self.x_extent(distance), self.y_extent(distance))

    # -----------------------------------------------------------------

    def pixel_area(self, distance):

        """
        This function ...
        :param distance: 
        :return: 
        """

        return self.x_extent(distance) * self.y_extent(distance)

    # -----------------------------------------------------------------

    def to_physical(self, distance, unit="pc"):

        """
        Thisn function ...
        :param distance:
        :param unit:
        :return:
        """

        return PhysicalPixelscale.from_angular(self, distance, unit=unit)

# -----------------------------------------------------------------

def only_angle(quantity, unit=None):

    """
    This function ...
    :param quantity:
    :param unit:
    :return:
    """

    if not hasattr(quantity, "unit"): return Angle(quantity, unit=unit)

    if "pix" in quantity.unit.bases: return quantity * Unit("pix")
    elif len(quantity.unit.bases) == 1: return quantity
    else: raise ValueError("Don't know what to do with " + str(quantity) + " to convert to angular pixelscale")

# -----------------------------------------------------------------

class PhysicalPixelscale(PhysicalExtent):

    """
    THis class ...
    """

    def __init__(self, x, y=None, unit=None):

        """
        The constructor ...
        """

        # Set y
        if y is None: y = x

        # Convert to just arcsec or just
        x = only_physical(x, unit=unit)
        y = only_physical(y, unit=unit)

        # Call the constructor of the base class
        super(PhysicalPixelscale, self).__init__(x, y)

    # -----------------------------------------------------------------

    @classmethod
    def from_angular(cls, pixelscale, distance, unit="pc"):

        """
        This function ...
        :param pixelscale:
        :param distance:
        :param unit:
        :return:
        """

        return cls(pixelscale.x_extent(distance, unit=unit), pixelscale.y_extent(distance, unit=unit))

    # -----------------------------------------------------------------

    @classmethod
    def from_shape_and_field(cls, shape, field):

        """
        This function ...
        :param shape: PixelShape
        :param field:
        :return:
        """

        x = field.x / float(shape.x)
        y = field.y / float(shape.y)
        return cls(x=x, y=y)

    # -----------------------------------------------------------------

    @property
    def abs_x(self):

        """
        This function ...
        :return:
        """

        return abs(self.x)

    # -----------------------------------------------------------------

    @property
    def abs_y(self):

        """
        This function ...
        :return:
        """

        return abs(self.y)

    # -----------------------------------------------------------------

    @property
    def average(self):

        """
        This function ...
        :return:
        """

        x_pixelscale = abs(self.x.to("pc"))
        y_pixelscale = abs(self.y.to("pc"))

        # Check if the x pixelscale and y pixelscale are close
        if not np.isclose(x_pixelscale.value, y_pixelscale.value, rtol=0.001):

            from ...core.tools.stringify import tostr

            # Warn about the difference in x and y pixelscale
            from ...core.basics.log import log
            log.warning("Averaging the pixelscale over the x and y direction may not be a good approximation:")
            log.warning("  * x pixelscale (absolute value) = " + tostr(x_pixelscale))
            log.warning("  * y pixelscale (absolute value) = " + tostr(y_pixelscale))

        # Return a single value for the pixelscale in parsec
        return 0.5 * (x_pixelscale + y_pixelscale)

    # -----------------------------------------------------------------

    def x_angle(self, distance, unit="arcsec"):

        """
        This function ...
        :param distance:
        :param unit:
        :return:
        """

        angular = (self.abs_x / distance).to(unit, equivalencies=dimensionless_angles())
        return angular

    # -----------------------------------------------------------------

    def y_angle(self, distance, unit="arcsec"):

        """
        This function ...
        :param distance:
        :param unit:
        :return:
        """

        angular = (self.abs_y / distance).to(unit, equivalencies=dimensionless_angles())
        return angular

    # -----------------------------------------------------------------

    def x_extent(self, unit="kpc"):

        """
        This function ...
        :param unit:
        :return:
        """

        return self.abs_x.to(unit)

    # -----------------------------------------------------------------

    def y_extent(self, unit="kpc"):

        """
        This function ...
        :param unit:
        :return:
        """

        return self.abs_y.to(unit)

    # -----------------------------------------------------------------

    @property
    def pixel_area(self):

        """
        This function ...
        :return:
        """

        return self.x_extent() * self.y_extent()

    # -----------------------------------------------------------------

    def to_angular(self, distance, unit="arcsec"):

        """
        Thisn function ...
        :param distance:
        :param unit:
        :return:
        """

        return Pixelscale.from_physical(self, distance, unit=unit)

# -----------------------------------------------------------------

def only_physical(quantity, unit=None):

    """
    This function ...
    :param quantity:
    :param unit:
    :return:
    """

    if not hasattr(quantity, "unit"): return quantity * Unit(unit)

    if "pix" in quantity.unit.bases: return quantity * Unit("pix")
    elif len(quantity.unit.bases) == 1: return quantity
    else: raise ValueError("Don't know what to do with " + str(quantity) + " to convert to physical pixelscale")

# -----------------------------------------------------------------
