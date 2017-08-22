#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.region.polygon Contains the PolygonRegion class and subclasses.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.coordinates import frame_transform_graph, SkyCoord

# Import the relevant PTS classes and modules
from .region import Region, PixelRegion, SkyRegion, PhysicalRegion
from ..basics.coordinate import PixelCoordinate, SkyCoordinate, PhysicalCoordinate
from .region import add_info, make_polygon_template, coordsys_name_mapping

# -----------------------------------------------------------------

class PolygonRegion(Region):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Set the points
        self.points = list(args)

        # Call the constructor of the base class
        super(PolygonRegion, self).__init__(**kwargs)

# -----------------------------------------------------------------

class PixelPolygonRegion(PolygonRegion, PixelRegion):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Verify the points
        for arg in args:
            if not isinstance(arg, PixelCoordinate): raise ValueError("Points should be of type 'PixelCoordinate'")

        # Call the constructor of the base class
        super(PixelPolygonRegion, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def add_point(self, point):

        """
        This function ...
        :return:
        """

        if not isinstance(point, PixelCoordinate): raise ValueError("Point should be of type 'PixelCoordinate'")

        # Add the point
        self.points.append(point)

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

        v = self.points
        coords = [(x, y) for x, y in zip(v.x, v.y)]
        val = "{:" + fmt + "}"
        temp = [val.format(x) for _ in coords for x in _]
        c = ",".join(temp)

        string = prefix + make_polygon_template().format(**locals())
        string = add_info(string, self)
        return string

    # -----------------------------------------------------------------

    def to_sky(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        return SkyPolygonRegion.from_pixel(self, wcs)

    # -----------------------------------------------------------------

    @classmethod
    def from_sky(cls, region, wcs):

        """
        This function ...
        :param region:
        :param wcs:
        :return:
        """

        points = []

        for coordinate in region.points:

            # Convert
            point = PixelCoordinate.from_sky(coordinate, wcs)

            # Add the point
            points.append(point)

        # Create the pixel circle region
        return cls(*points, meta=region.meta, label=region.label, include=region.include,
                   appearance=region.appearance)

# -----------------------------------------------------------------

class SkyPolygonRegion(PolygonRegion, SkyRegion):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Verify the points
        for index, arg in enumerate(args):
            if not isinstance(arg, SkyCoordinate): raise ValueError("Points should be of type 'SkyCoordinate'")

        # Call the constructor of the base class
        super(SkyPolygonRegion, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def add_point(self, point):

        """
        This function ...
        :param point:
        :return:
        """

        if not isinstance(point, SkyCoordinate): raise ValueError("Point should be of type 'SkyCoordinate'")

        # Add the point
        self.points.append(point)

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

        v = self.points.transform_to(frame)
        coords = [(x.to('deg').value, y.to('deg').value) for x, y in
                  zip(v.spherical.lon, v.spherical.lat)]
        val = "{:" + fmt + "}"
        temp = [val.format(x) for _ in coords for x in _]
        c = ",".join(temp)

        string = prefix + make_polygon_template().format(**locals())
        string = add_info(string, self)
        return string

    # -----------------------------------------------------------------

    @classmethod
    def from_pixel(cls, region, wcs):

        """
        This function ...
        :parma region:
        :param wcs:
        :return:
        """

        coordinates = []

        for point in region.points:
            coordinate = SkyCoordinate.from_pixel(point, wcs)
            coordinates.append(coordinate)

        # Create and return
        return cls(*coordinates, meta=region.meta, label=region.label, include=region.include,
                   appearance=region.appearance)

    # -----------------------------------------------------------------

    def to_pixel(self, wcs):

        """
        Thisn function ...
        :param wcs:
        :return:
        """

        return PixelPolygonRegion.from_sky(self, wcs)

# -----------------------------------------------------------------

class PhysicalPolygonRegion(PolygonRegion, PhysicalRegion):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Verify the points
        for arg in args:
            if not isinstance(arg, PhysicalCoordinate): raise ValueError("Points should be of type 'PhysicalCoordinate'")

        # Call the constructor of the base class
        super(PhysicalPolygonRegion, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def add_point(self, point):

        """
        This function ...
        :param point:
        :return:
        """

        if not isinstance(point, PhysicalCoordinate): raise ValueError("Point should be of type 'PhysicalCoordinate'")

        # Add the point
        self.points.append(point)

# -----------------------------------------------------------------
