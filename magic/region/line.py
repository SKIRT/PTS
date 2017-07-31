#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.region.line Contains the LineRegion class and subclasses.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.units import Unit
from astropy.coordinates import frame_transform_graph

# Import the relevant PTS classes and modules
from .region import Region, PixelRegion, SkyRegion, PhysicalRegion
from ..basics.coordinate import PixelCoordinate, SkyCoordinate, PhysicalCoordinate
from ..basics.stretch import PixelStretch, SkyStretch, PhysicalStretch
from .rectangle import PixelRectangleRegion, SkyRectangleRegion, PhysicalRectangleRegion
from ..tools import coordinates
from .region import make_line_template, add_info, coordsys_name_mapping

# -----------------------------------------------------------------

class LineRegion(Region):

    """
    This class ...
    """

    def __init__(self, start, end, **kwargs):

        """
        The constructor ...
        :param start:
        :param end:
        :param kwargs:
        """

        # Set the attributes
        self.start = start
        self.end = end

        # Call the constructor of the base class
        super(LineRegion, self).__init__(**kwargs)

    # -----------------------------------------------------------------

    @property
    def axis1_center(self):

        """
        This property ...
        :return:
        """

        return 0.5*(self.start.axis1 + self.end.axis1)

    # -----------------------------------------------------------------

    @property
    def axis2_center(self):

        """
        This property ...
        :return:
        """

        return 0.5*(self.start.axis2 + self.end.axis2)

    # -----------------------------------------------------------------

    @property
    def axis1_min(self):

        """
        This function ...
        :return:
        """

        return min(self.start.axis1, self.end.axis1)

    # -----------------------------------------------------------------

    @property
    def axis1_max(self):

        """
        This function ...
        :return:
        """

        return max(self.start.axis1, self.end.axis1)

    # -----------------------------------------------------------------

    @property
    def axis2_min(self):

        """
        This function ...
        :return:
        """

        return min(self.start.axis2, self.end.axis2)

    # -----------------------------------------------------------------

    @property
    def axis2_max(self):

        """
        This function ...
        :return:
        """

        return max(self.start.axis2, self.end.axis2)

# -----------------------------------------------------------------

class PixelLineRegion(LineRegion, PixelRegion):

    """
    This class ...
    """

    def __init__(self, start, end, **kwargs):

        """
        The constructor ...
        :param start:
        :param end:
        :param kwargs:
        """

        # Check the arguments
        if not isinstance(start, PixelCoordinate): raise ValueError("Start must be a pixel coordinate")
        if not isinstance(end, PixelCoordinate): raise ValueError("End must be a pixel coordinate")

        # Call the constructor of the base class
        super(PixelLineRegion, self).__init__(start, end, **kwargs)

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
        This setter ...
        :param value:
        :return:
        """

        offset = value - self.center

        self.start += offset
        self.end += offset

    # -----------------------------------------------------------------

    @property
    def radius(self):

        """
        This function ...
        :return:
        """

        x_min = self.axis1_min
        x_max = self.axis1_max
        y_min = self.axis2_min
        y_max = self.axis2_max

        return PixelStretch(0.5 * (x_max - x_min), 0.5 * (y_max - y_min))

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        if self.include: prefix = ""
        else: prefix = "-"

        fmt = '.4f'

        x1 = self.start.x
        y1 = self.end.y

        x2 = self.end.x
        y2 = self.end.y

        string = prefix + make_line_template(fmt).format(**locals())
        string = add_info(string, self)
        return string

# -----------------------------------------------------------------

class SkyLineRegion(LineRegion, SkyRegion):

    """
    This class ...
    """

    def __init__(self, start, end, **kwargs):

        """
        The constructor ...
        :param start:
        :param end:
        :param kwargs:
        """

        # Check the arguments
        if not isinstance(start, SkyCoordinate): raise ValueError("Start must be a sky coordinate")
        if not isinstance(end, SkyCoordinate): raise ValueError("End must be a sky coordinate")

        # Call the constructor of the base class
        super(SkyLineRegion, self).__init__(start, end, **kwargs)

    # -----------------------------------------------------------------

    @property
    def center(self):

        """
        This function ...
        :return:
        """

        return SkyCoordinate(self.axis1_center, self.axis2_center)

    # -----------------------------------------------------------------

    @property
    def radius(self):

        """
        This function ...
        :return:
        """

        dec_start = self.start.dec.to("deg").value
        dec_end = self.end.dec.to("deg").value

        dec_center = 0.5 * (dec_start + dec_end)

        ra_start = self.start.ra.to("deg").value
        ra_end = self.end.ra.to("deg").value

        # Calculate the actual RA and DEC distance in degrees
        ra_distance = abs(coordinates.ra_distance(dec_center, ra_start, ra_end))
        dec_distance = abs(dec_end - dec_start)

        ra_span = ra_distance * Unit("deg")
        dec_span = dec_distance * Unit("deg")

        return SkyStretch(0.5 * ra_span, 0.5 * dec_span)

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        coordsys = "fk5"
        # convert coordsys string to coordsys object
        if coordsys in coordsys_name_mapping: frame = frame_transform_graph.lookup_name(coordsys_name_mapping[coordsys])
        else: frame = None  # for pixel/image/physical frames

        fmt = '.4f'

        if self.include: prefix = ""
        else: prefix = "-"

        x1 = float(self.start.transform_to(frame).spherical.lon.to('deg').value)
        y1 = float(self.start.transform_to(frame).spherical.lat.to('deg').value)

        x2 = float(self.end.transform_to(frame).spherical.lon.to('deg').value)
        y2 = float(self.end.transform_to(frame).spherical.lat.to('deg').value)

        # TO HAVE IT IN hh:mm:ss NOTATION:

        # print(vars(reg.start.transform_to(frame)))
        # skycoordinate = SkyCoordinate(reg.start.transform_to(frame).spherical.lon, reg.start.transform_to(frame).spherical.lat, frame=frame, representation="spherical")
        # str1 = skycoordinate.to_string('hmsdms').replace("d", ":").replace("h", ":").replace("m", ":").replace("s ", ",")[:-1]

        # skycoordinate = SkyCoordinate(reg.end.transform_to(frame).spherical.lon, reg.end.transform_to(frame).spherical.lat, frame=frame, representation="spherical")
        # str2 = skycoordinate.to_string('hmsdms').replace("d", ":").replace("h", ":").replace("m", ":").replace("s ", ",")[:-1]

        string = prefix + make_line_template(fmt).format(**locals())
        string = add_info(string, self)
        return string

# -----------------------------------------------------------------

class PhysicalLineRegion(LineRegion, PhysicalRegion):

    """
    This class ...
    """

    def __init__(self, start, end, **kwargs):

        """
        The constructor ...
        :param start:
        :param end:
        :param kwargs:
        """

        # Check the arguments
        if not isinstance(start, PhysicalCoordinate): raise ValueError("Start must be a physical coordinate")
        if not isinstance(end, PhysicalCoordinate): raise ValueError("End must be a physical coordinate")

        # Call the constructor of the base class
        super(PhysicalLineRegion, self).__init__(start, end, **kwargs)

    # -----------------------------------------------------------------

    @property
    def center(self):

        """
        This function ...
        :return:
        """

        return PhysicalCoordinate(self.axis1_center, self.axis2_center)

# -----------------------------------------------------------------
