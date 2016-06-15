#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.basics.skyregion Contains the SkyRegion class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import astronomical modules
import pyregion
from astropy.units import Unit
from astropy.coordinates import Angle

# Import the relevant PTS classes and modules
from .vector import Extent
from .skygeometry import SkyCoordinate, SkyLine, SkyCircle, SkyEllipse, SkyRectangle, SkyPolygon

# -----------------------------------------------------------------

class SkyRegion(list):

    """
    This class ...
    """

    def __init__(self):

        """
        This function ...
        :return:
        """

        # Call the constructor of the base class
        super(SkyRegion, self).__init__()

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, only=None, ignore=None, color=None, ignore_color=None):

        """
        This function ...
        :param path:
        :param only:
        :param ignore:
        :param color:
        :param ignore_color:
        :return:
        """

        # Create a new sky region
        region = cls()

        # Open the region file with pyregion and check if its in sky coordinates
        _region = pyregion.open(path)
        def check_sky_coord(_region):
            if [s for s in _region if s.coord_format != "fk5"]: return False
            else: return True
        if not check_sky_coord(_region): raise ValueError("Region is not in sky coordinates")

        # Loop over all shapes in the region
        for shape in _region:

            # Meta information
            meta = {}
            if "text" in shape.attr[1]: meta["text"] = shape.attr[1]["text"]
            if "color" in shape.attr[1]: meta["color"] = shape.attr[1]["color"]
            if "point" in shape.attr[1]: meta["point"] = shape.attr[1]["point"]

            # Check the color of the shape
            if color is not None and shape.attr[1]["color"] != color: continue
            if ignore_color is not None and shape.attr[1]["color"] == ignore_color: continue

            # Get the coordinate format of this shape
            coord_format = shape.coord_format

            # The shape is a point -> SkyCoord
            if shape.name == "point":

                # Get RA and declination
                ra = shape.coord_list[0]
                dec = shape.coord_list[1]

                # Create sky coordinate
                coordinate = SkyCoordinate(ra=ra, dec=dec, unit="deg", frame=coord_format, meta=meta)

                new_shape = coordinate

            # The shape is a line or vector -> SkyLine
            elif shape.name == "line" or shape.name == "vector":

                # Get the RA and declination of the two points
                ra_1 = shape.coord_list[0]
                dec_1 = shape.coord_list[1]
                ra_2 = shape.coord_list[2]
                dec_2 = shape.coord_list[3]

                # Create the sky coordinates
                coordinate_1 = SkyCoordinate(ra=ra_1, dec=dec_1, unit="deg", frame=coord_format)
                coordinate_2 = SkyCoordinate(ra=ra_2, dec=dec_2, unit="deg", frame=coord_format)

                # Create the SkyLine object
                line = SkyLine(coordinate_1, coordinate_2, meta=meta)

                new_shape = line

            # The shape is a circle -> SkyCircle
            elif shape.name == "circle":

                # Get the RA and declination of the center and the radius
                ra_center = shape.coord_list[0]
                dec_center = shape.coord_list[1]
                radius = shape.coord_list[2] * Unit("deg")

                # Create a sky cooridnate for the center
                center = SkyCoordinate(ra=ra_center, dec=dec_center, unit="deg", frame=coord_format)

                # Create a SkyCircle object and add it to the region
                circle = SkyCircle(center, radius, meta=meta)

                new_shape = circle

            # The shape is an ellipse -> SkyEllipse
            elif shape.name == "ellipse":

                # Get the RA and declination of the center
                ra_center = shape.coord_list[0]
                dec_center = shape.coord_list[1]
                center = SkyCoordinate(ra=ra_center, dec=dec_center, unit="deg", frame=coord_format)

                # Get the radius
                x_radius = shape.coord_list[2] * Unit("deg")
                y_radius = shape.coord_list[3] * Unit("deg")
                radius = Extent(x_radius, y_radius)

                # Get the angle
                angle = Angle(shape.coord_list[4], "deg")

                # Create a SkyEllipse object and add it to the region
                ellipse = SkyEllipse(center, radius, angle, meta=meta)

                new_shape = ellipse

            # The shape is a rectangle -> SkyRectangle
            elif shape.name == "box":

                # Get the RA and declination of the center
                ra_center = shape.coord_list[0]
                dec_center = shape.coord_list[1]
                center = SkyCoordinate(ra=ra_center, dec=dec_center, unit="deg", frame=coord_format)

                # Get the width and height
                width = shape.coord_list[2] * Unit("deg")
                height = shape.coord_list[3] * Unit("deg")

                # Create radius
                radius = Extent(0.5 * width, 0.5*height)

                # Get the angle
                angle = Angle(shape.coord_list[4], "deg")

                # Create a SkyRectangle and add it to the region
                rectangle = SkyRectangle(center, radius, angle, meta=meta)

                new_shape = rectangle

            # The shape is a polygon -> SkyPolygon
            elif shape.name == "polygon":

                number_of_points = 0.5 * len(shape.coord_list)
                assert int(number_of_points) == number_of_points
                number_of_points = int(number_of_points)

                # Create a new SkyPolygon
                polygon = SkyPolygon(meta=meta)

                # Get the RA and declination of the different points
                for i in range(number_of_points):

                    # Create a new SkyCoordinate object
                    ra = shape.coord_list[2*i]
                    dec = shape.coord_list[2*i + 1]
                    coordinate = SkyCoordinate(ra=ra, dec=dec, unit="deg", frame=coord_format)

                    # Add the coordinate to the polygon
                    polygon.add_point(coordinate)

                new_shape = polygon

            # Unrecognized shape
            else: raise ValueError("Unrecognized shape (should be point, line, vector, circle, ellipse, box or polygon")

            region.append(new_shape)

        # Return the new region
        return region

    # -----------------------------------------------------------------

    def append(self, shape):

        """
        This function ...
        :param shape:
        :return:
        """

        # Check whether the shape is in sky coordinates
        if not shape.__class__.__name__.startswith("Sky"): raise ValueError("Shape must be SkyCoordinate, SkyLine, SkyCircle, SkyEllipse, SkyRectangle, SkyPolygon or SkyComposite")

        # Otherwise, add the shape
        super(SkyRegion, self).append(shape)

    # -----------------------------------------------------------------

    def __mul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        # Create a new region
        new_region = SkyRegion()
        for shape in self: new_region.append(shape * value)

        # Return the new region
        return new_region

    # -----------------------------------------------------------------

    def __truediv__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self.__div__(value)

    # -----------------------------------------------------------------

    def __div__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        # Create a new region
        new_region = SkyRegion()
        for ellipse in self: new_region.append(ellipse / value)

        # Return the new region
        return new_region

    # -----------------------------------------------------------------

    def to_pixel(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        # Avoid circular import at module level
        from .region import Region

        # Create a new region
        region = Region()

        # Fill the new list
        for shape in self: region.append(shape.to_pixel(wcs))

        # Return region in pixel coordinates
        return region

    # -----------------------------------------------------------------

    def save(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Create a file
        f = open(path, 'w')

        # Initialize the region string
        print("# Region file format: DS9 version 4.1", file=f)

        # Write the coordinate system
        print("fk5\n", file=f)

        # Loop over all shapes, get string and print it to the region file
        for shape in self: print(shape.to_region_string(coordinate_system=False), file=f)

        # Close the file
        f.close()

# -----------------------------------------------------------------
