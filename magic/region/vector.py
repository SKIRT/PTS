#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.region.text Contains the TextRegion class and subclasses.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.coordinates import Angle
from astropy.units import Quantity
from astropy.coordinates import frame_transform_graph

# Import the relevant PTS classes and modules
from .region import Region, PixelRegion, SkyRegion, PhysicalRegion
from ..basics.coordinate import PixelCoordinate, SkyCoordinate, PhysicalCoordinate
from ..basics.stretch import PixelStretch, SkyStretch, PhysicalStretch
from .rectangle import PixelRectangleRegion, SkyRectangleRegion
from ..basics.mask import Mask
from .region import make_vector_template, add_info, coordsys_name_mapping

# -----------------------------------------------------------------

class VectorRegion(Region):

    """
    This class ...
    """

    def __init__(self, start, length, angle, **kwargs):

        """
        The constructor ...
        :param start:
        :param length:
        :param angle:
        :param kwargs:
        """

        # Check the angle
        if not isinstance(angle, Angle): raise ValueError("Angle must be a Astropy Angle object")

        # Set the attributes
        self.start = start
        self.length = length
        self.angle = angle

        # Call the constructor of the base class
        super(VectorRegion, self).__init__(**kwargs)

    # -----------------------------------------------------------------

    @property
    def axis1_center(self):

        """
        This property ...
        :return:
        """

        raise NotImplementedError("Not implemented yet")

    # -----------------------------------------------------------------

    @property
    def axis2_center(self):

        """
        This property ...
        :return:
        """

        raise NotImplementedError("Not implemented yet")

# -----------------------------------------------------------------

class PixelVectorRegion(VectorRegion, PixelRegion):

    """
    This class ...
    """

    def __init__(self, start, length, angle, **kwargs):

        """
        This function ...
        :param start:
        :param length:
        :param angle:
        :param kwargs:
        """

        # Check the start coordinate
        if not isinstance(start, PixelCoordinate): raise ValueError("Start must be pixel coordinate")

        # Check the length
        if not isinstance(length, float): raise ValueError("Length must be float")

        # Call the constructor of VectorRegion class
        VectorRegion.__init__(self, start, length, angle, **kwargs)

    # -----------------------------------------------------------------

    @property
    def center(self):

        """
        This function ...
        :return:
        """

        return PixelCoordinate(self.axis1_center, self.axis2_center)

    # -----------------------------------------------------------------

    @center.setter
    def center(self, value):

        """
        This function ...
        :return:
        """

        offset = value - self.center

        self.start += offset

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        if self.include: prefix = ""
        else: prefix = "-"

        coordsys = 'fk5'
        fmt = '.4f'
        radunit = 'deg'

        if radunit == 'arcsec':
            if coordsys in coordsys_name_mapping.keys(): radunitstr = '"'
            else: raise ValueError('Radius unit arcsec not valid for coordsys {}'.format(coordsys))
        else: radunitstr = ''

        x = self.start.x
        y = self.start.y
        l = self.length
        ang = self.angle.to("deg").value

        string = prefix + make_vector_template(fmt, radunitstr).format(**locals())
        string = add_info(string, self)
        return string

# -----------------------------------------------------------------

class SkyVectorRegion(VectorRegion, SkyRegion):

    """
    This class ...
    """

    def __init__(self, start, length, angle, **kwargs):

        """
        This function ...
        :param start:
        :param length:
        :param angle:
        :param kwargs:
        """

        # Check the start coordinate
        if not isinstance(start, SkyCoordinate): raise ValueError("Start must be sky coordinate")

        # Check the length
        if not isinstance(length, Quantity): raise ValueError("Length must be an angular quantity")

        # Call the constructor of VectorRegion class
        VectorRegion.__init__(self, start, length, angle, **kwargs)

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        if self.include: prefix = ""
        else: prefix = "-"

        coordsys = 'fk5'
        fmt = '.4f'
        radunit = 'deg'

        if radunit == 'arcsec':
            if coordsys in coordsys_name_mapping.keys(): radunitstr = '"'
            else: raise ValueError('Radius unit arcsec not valid for coordsys {}'.format(coordsys))
        else: radunitstr = ''

        # convert coordsys string to coordsys object
        if coordsys in coordsys_name_mapping: frame = frame_transform_graph.lookup_name(coordsys_name_mapping[coordsys])
        else: frame = None  # for pixel/image/physical frames

        x = float(self.start.transform_to(frame).spherical.lon.to('deg').value)
        y = float(self.start.transform_to(frame).spherical.lat.to('deg').value)
        l = float(self.length.to(radunit).value)
        ang = float(self.angle.to('deg').value)

        string = prefix + make_vector_template(fmt, radunitstr).format(**locals())
        string = add_info(string, self)
        return string

# -----------------------------------------------------------------

class PhysicalVectorRegion(VectorRegion, PhysicalRegion):

    """
    This class ...
    """

    def __init__(self, start, length, angle, **kwargs):

        """
        This function ...
        :param start:
        :param length:
        :param angle:
        :param kwargs:
        """

        # Check the start coordinate
        if not isinstance(start, PhysicalCoordinate): raise ValueError("Start must be physical coordinate")

        # Check the length
        if not isinstance(length, Quantity): raise ValueError("Length must be a physical quantity of length")

        # Call the constructor of VectorRegion class
        VectorRegion.__init__(self, start, length, angle, **kwargs)

# -----------------------------------------------------------------
