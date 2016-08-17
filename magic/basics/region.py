#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.basics.region Contains the Region class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import astronomical modules
import pyregion
from astropy.coordinates import Angle

# Import the relevant PTS classes and modules
from .vector import Position, Extent
from .geometry import Coordinate, Line, Circle, Ellipse, Rectangle, Polygon, Composite
from .mask import Mask

# -----------------------------------------------------------------

class Region(list):

    """
    This class ...
    """

    default_extension = "reg"

    # -----------------------------------------------------------------

    def __init__(self):

        """
        This function ...
        :return:
        """

        # Call the constructor of the base class
        super(Region, self).__init__()

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

        # Create a new region
        region = cls()

        # Open the region file with pyregion and check if its in image coordinates
        try:
            _region = pyregion.open(path)
            if not _region.check_imagecoord(): raise IOError("Region is not in image coordinates")
        except ValueError: # If a ValueError comes out, assume the region file is empty (no shapes)
            _region = []

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

            # If the shape is a point -> Position
            if shape.name == "point":

                if only is not None and "point" not in only: continue
                if ignore is not None and "point" in ignore: continue

                # Get the position
                x = shape.coord_list[0]
                y = shape.coord_list[1]
                coordinate = Coordinate(x, y, meta=meta)

                # Add the position to the region
                new_shape = coordinate

            # If the shape is a line -> Line
            elif shape.name == "line" or shape.name == "vector":

                if only is not None and "line" not in only: continue
                if ignore is not None and "line" in ignore: continue

                # Get the position of the two points
                x_1 = shape.coord_list[0]
                y_1 = shape.coord_list[1]
                x_2 = shape.coord_list[2]
                y_2 = shape.coord_list[3]

                # Create the positions
                position_1 = Coordinate(x_1, y_1)
                position_2 = Coordinate(x_2, y_2)

                # Create the line
                line = Line(position_1, position_2, meta=meta)
                new_shape = line

            # If the shape is a circle -> Circle
            elif shape.name == "circle":

                if only is not None and "circle" not in only: continue
                if ignore is not None and "circle" in ignore: continue

                # Get the position of the center
                x_center = shape.coord_list[0]
                y_center = shape.coord_list[1]
                center = Coordinate(x_center, y_center)

                # Get the radius
                radius = shape.coord_list[2]

                # Create a circle
                circle = Circle(center, radius, meta=meta)
                new_shape = circle

            # If the shape is an ellipse -> Ellipse
            elif shape.name == "ellipse":

                if only is not None and "ellipse" not in only: continue
                if ignore is not None and "ellipse" in ignore: continue

                # Get the position of the center
                x_center = shape.coord_list[0]
                y_center = shape.coord_list[1]
                center = Coordinate(x_center, y_center)

                # Get the radius
                x_radius = shape.coord_list[2]
                y_radius = shape.coord_list[3]
                radius = Extent(x_radius, y_radius)

                # Get the angle
                angle = Angle(shape.coord_list[4], "deg")

                # Create an ellipse
                ellipse = Ellipse(center, radius, angle, meta=meta)
                new_shape = ellipse

            # If the shape is a rectangle -> Rectangle
            elif shape.name == "box":

                if only is not None and "box" not in only: continue
                if ignore is not None and "box" in ignore: continue

                # Get the position of the center
                x_center = shape.coord_list[0]
                y_center = shape.coord_list[1]
                center = Coordinate(x_center, y_center)

                # Get the width and height
                width = shape.coord_list[2]
                height = shape.coord_list[3]

                # Create radius
                radius = Extent(0.5 * width, 0.5 * height)

                # Get the angle
                angle = Angle(shape.coord_list[4], "deg")

                # Create a Rectangle
                rectangle = Rectangle(center, radius, angle, meta=meta)
                new_shape = rectangle

            # If the shape is a polygon -> Polygon
            elif shape.name == "polygon":

                if only is not None and "polygon" not in only: continue
                if ignore is not None and "polygon" in ignore: continue

                # Get the number of points in the polygon
                number_of_points = 0.5 * len(shape.coord_list)
                assert int(number_of_points) == number_of_points
                number_of_points = int(number_of_points)

                # Create a new Polygon
                polygon = Polygon(meta=meta)

                # Get the position of the different points
                for i in range(number_of_points):

                    # Create a new Position
                    x = shape.coord_list[2*i]
                    y = shape.coord_list[2*i + 1]
                    position = Coordinate(x, y)

                    # Add the point to the polygon
                    polygon.add_point(position)

                # Add the polygon to the region
                new_shape = polygon

            # Unrecognized shape
            else: raise ValueError("Unrecognized shape: " + shape.name)

            # If this shape should be excluded
            if shape.exclude:

                previous_shape = region[len(region)-1]
                region[len(region)-1] = Composite(previous_shape, new_shape, meta=previous_shape.meta)

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

        # Check whether the argument is a valid shape
        if not (shape.__class__.__name__ == "Coordinate" or shape.__class__.__name__ == "Line"
                or shape.__class__.__name__ == "Circle" or shape.__class__.__name__ == "Ellipse"
                or shape.__class__.__name__ == "Rectangle" or shape.__class__.__name__ == "Polygon"
                or shape.__class__.__name__ == "Composite"):
            raise ValueError("Shape must of of type Coordinate, Line, Circle, Ellipse, Rectangle, Polygon or Composite")

        # Otherwise, add the shape
        super(Region, self).append(shape)

    # -----------------------------------------------------------------

    def rectangles(self):

        """
        This function ...
        :return:
        """

        return [shape for shape in self if isinstance(shape, Rectangle)]

    # -----------------------------------------------------------------

    def circles(self):

        """
        This function ...
        :return:
        """

        return [shape for shape in self if isinstance(shape, Circle)]

    # -----------------------------------------------------------------

    def ellipses(self):

        """
        This function ...
        :return:
        """

        return [shape for shape in self if isinstance(shape, Ellipse)]

    # -----------------------------------------------------------------

    def cropped(self, x_min, x_max, y_min, y_max):

        """
        This function ...
        :param x_min:
        :param x_max:
        :param y_min:
        :param y_max:
        :return:
        """

        # Create a new region
        region = Region()

        # Loop over all shapes
        for shape in self:

            # Ignore shapes that fall outside of the crop region
            if shape.x < x_min or shape.x > x_max or shape.y < y_min or shape.y > y_max: continue

            # Create a new shape that has its coordinates shifted
            new_shape = shape - Extent(x_min, y_min)

            # Add the new shape to the new region
            region.append(new_shape)

        # Return the new region
        return region

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

    def to_sky(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        # Avoid circular import at module scope
        from .skyregion import SkyRegion

        # Create a new SkyRegion
        region = SkyRegion()

        # Add the shapes to the sky region
        for shape in self: region.append(shape.to_sky(wcs))

        # Return the region
        return region

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

    def to_mpl_patches(self):

        """
        This function ...
        :return:
        """

        # Initialize a list to contain the patches
        patches = []

        # Add the patches from the shapes
        for shape in self: patches.append(shape.to_mpl_patch())

        # Return the list of patches
        return patches

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

        # Print the coordinate system
        print("image", file=f)

        # Loop over all shapes, get string and print it to the region file
        for shape in self: print(shape.to_region_string(coordinate_system=False), file=f)

        # Close the file
        f.close()

    # -----------------------------------------------------------------

    def homogenized(self):

        """
        This function returns a copy of the region where Composite objects have been dissolved into their individual components
        :return:
        """

        import copy

        new = Region()

        for shape in self:

            if isinstance(shape, Composite):

                copy_base = copy.deepcopy(shape.base)
                copy_exclude = copy.deepcopy(shape.exclude)
                copy_base.meta = shape.meta
                copy_exclude.meta = shape.meta

                new.append(copy_base)
                new.append(copy_exclude)

            else: new.append(copy.deepcopy(shape))

        # Return the new region
        return new

# -----------------------------------------------------------------
