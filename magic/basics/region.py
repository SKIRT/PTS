#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.basics.region Contains the Region class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
import pyregion
from astropy import units as u
from astropy.coordinates import Angle

# Import the relevant AstroMagic classes and modules
from .vector import Position, Extent
from .geometry import Line, Circle, Ellipse, Rectangle, Polygon, Composite
from .skygeometry import SkyCoord, SkyLine, SkyCircle, SkyEllipse, SkyRectangle, SkyPolygon
from .mask import Mask

# -----------------------------------------------------------------

class Region(list):

    """
    This class ...
    """

    def __init__(self):

        """
        This function ...
        :return:
        """

        # Call the constructor of the base class
        super(Region, self).__init__()

        # List of shapes to exclude
        self.exclude = []

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, only=None, ignore=None):

        """
        This function ...
        :param path:
        :param only:
        :param ignore:
        :return:
        """

        # Create a new region
        region = cls()

        # Open the region file with pyregion and check if its in image coordinates
        _region = pyregion.open(path)
        if not _region.check_imagecoord(): raise ValueError("Region is not in image coordinates")

        # Loop over all shapes in the region
        for shape in _region:

            # If the shape is a point -> Position
            if shape.name == "point":

                if only is not None and "point" not in only: continue

                # Get the position
                x = shape.coord_list[0]
                y = shape.coord_list[1]
                position = Position(x, y)

                # Add the position to the region
                #region.append(position)
                new_shape = position

            # If the shape is a line -> Line
            elif shape.name == "line" or shape.name == "vector":

                if only is not None and "line" not in only and "vector" not in only: continue

                # Get the position of the two points
                x_1 = shape.coord_list[0]
                y_1 = shape.coord_list[1]
                x_2 = shape.coord_list[2]
                y_2 = shape.coord_list[3]

                # Create the positions
                position_1 = Position(x_1, y_1)
                position_2 = Position(x_2, y_2)

                # Create the line and add it to the region
                line = Line(position_1, position_2)
                #region.append(line)
                new_shape = line

            # If the shape is a circle -> Circle
            elif shape.name == "circle":

                if only is not None and "circle" not in only: continue

                # Get the position of the center
                x_center = shape.coord_list[0]
                y_center = shape.coord_list[1]
                center = Position(x_center, y_center)

                # Get the radius
                radius = shape.coord_list[2]

                # Create a circle and add it to the region
                circle = Circle(center, radius)
                #region.append(circle)
                new_shape = circle

            # If the shape is an ellipse -> Ellipse
            elif shape.name == "ellipse":

                if only is not None and "ellipse" not in only: continue

                # Get the position of the center
                x_center = shape.coord_list[0]
                y_center = shape.coord_list[1]
                center = Position(x_center, y_center)

                # Get the radius
                x_radius = shape.coord_list[2]
                y_radius = shape.coord_list[3]
                radius = Extent(x_radius, y_radius)

                # Get the angle
                angle = Angle(shape.coord_list[4], "deg")

                # Create an ellipse and add it to the region
                ellipse = Ellipse(center, radius, angle)
                #region.append(ellipse)
                new_shape = ellipse

            # If the shape is a rectangle -> Rectangle
            elif shape.name == "box":

                if only is not None and "box" not in only: continue

                # Get the position of the center
                x_center = shape.coord_list[0]
                y_center = shape.coord_list[1]
                center = Position(x_center, y_center)

                # Get the width and height
                width = shape.coord_list[2]
                height = shape.coord_list[3]

                # Create radius
                radius = Extent(0.5 * width, 0.5 * height)

                # Create a Rectangle and add it to the region
                rectangle = Rectangle(center, radius)
                #region.append(rectangle)
                new_shape = rectangle

            # If the shape is a polygon -> Polygon
            elif shape.name == "polygon":

                if only is not None and "polygon" not in only: continue

                # Get the number of points in the polygon
                number_of_points = 0.5 * len(shape.coord_list)
                assert int(number_of_points) == number_of_points
                number_of_points = int(number_of_points)

                # Create a new Polygon
                polygon = Polygon()

                # Get the position of the different points
                for i in range(number_of_points):

                    # Create a new Position
                    x = shape.coord_list[2*i]
                    y = shape.coord_list[2*i + 1]
                    position = Position(x, y)

                    # Add the point to the polygon
                    polygon.add_point(position)

                # Add the polygon to the region
                #region.append(polygon)
                new_shape = polygon

            # Unrecognized shape
            else: raise ValueError("Unrecognized shape")

            # If this shape should be excluded
            if shape.exclude:

                previous_shape = region[len(region)-1]
                region[len(region)-1] = Composite(previous_shape, new_shape)

            # Add the shape to the region
            else: region.append(new_shape)

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
        if not (shape.__class__.__name__ == "Position" or shape.__class__.__name__ == "Line"
                or shape.__class__.__name__ == "Circle" or shape.__class__.__name__ == "Ellipse"
                or shape.__class__.__name__ == "Rectangle" or shape.__class__.__name__ == "Polygon"):
            raise ValueError("Shape must of of type Position, Line, Circle, Ellipse, Rectangle or Polygon")

        # Otherwise, add the shape
        super(Region, self).append(shape)

    # -----------------------------------------------------------------

    def __mul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        new_region = Region()
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

        return self.__div__(value)

    # -----------------------------------------------------------------

    def __div__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        new_region = Region()
        for shape in self: new_region.append(shape / value)

        # Return the new region
        return new_region

    # -----------------------------------------------------------------

    def to_sky_coordinates(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        # Avoid circular import at module scope
        from .skyregion import SkyRegion

        # Initialize a new list contain the shapes in sky coordinates
        new_region = SkyRegion()

        # Fill the new list
        for shape in self:

            if shape.__class__.__name__ == "Position": new_region.append(SkyCoord.from_pixel(shape.x, shape.y, wcs, mode="wcs"))
            elif shape.__class__.__name__ == "Line": new_region.append(SkyLine.from_line(shape, wcs))
            elif shape.__class__.__name__ == "Circle": new_region.append(SkyCircle.from_circle(shape, wcs))
            elif shape.__class__.__name__ == "Ellipse": new_region.append(SkyEllipse.from_ellipse(shape, wcs))
            elif shape.__class__.__name__ == "Rectangle": new_region.append(SkyRectangle.from_rectangle(shape, wcs))
            elif shape.__class__.__name__ == "Polygon": new_region.append(SkyPolygon.from_polygon(shape, wcs))
            else: raise ValueError("Unrecognized shape")

        # Return the list of ellipses in image coordinates
        return new_region

    # -----------------------------------------------------------------

    def to_mask(self, x_size, y_size):

        """
        This function ...
        :param x_size:
        :param y_size:
        :return:
        """

        # Create empty mask
        mask = Mask.empty(x_size, y_size)

        # Add each shape to the mask
        for shape in self: mask += Mask.from_shape(shape, x_size, y_size)

        # Return the mask
        return mask

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

        # Loop over all ellipses
        for ellipse in self:

            # Get aperture properties
            center = ellipse.center
            major = ellipse.major
            minor = ellipse.minor
            angle = ellipse.angle.degree

            # Write to region file
            print("image;ellipse({},{},{},{},{})".format(center.x+1, center.y+1, major, minor, angle), file=f)

        # Close the file
        f.close()

# -----------------------------------------------------------------
