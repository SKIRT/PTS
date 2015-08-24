#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import standard modules
import copy

# Import astronomical modules
import pyregion

# *****************************************************************

class Region(object):

    """
    This class ...
    """

    # *****************************************************************

    def __init__(self, region):

        """
        The constructor ...
        :param region:
        :return:
        """

        # Set the internal pyregion object
        self._region = region

        # Set as unactive initially
        self.selected = False

    # *****************************************************************

    def select(self):

        """
        This function ...
        :return:
        """

        self.selected = True

    # *****************************************************************

    def deselect(self):

        """
        This function ...
        :return:
        """

        self.selected = False

    # *****************************************************************

    @property
    def number_of_shapes(self):

        """
        This function returns the number of shapes in this region
        :return:
        """

        return len(self._region)

# *****************************************************************

def ellipse_parameters(shape):

    """
    This function ...
    :param shape:
    :return:
    """

    x_center = shape.coord_list[0]
    y_center = shape.coord_list[1]
    x_radius = shape.coord_list[2]
    y_radius = shape.coord_list[3]
    
    return x_center, y_center, x_radius, y_radius

# *****************************************************************

def get_enclosing_box(shape):

    """
    This function ...
    :param shape:
    :return:
    """

    # TODO: make it work for shapes other than ellipses!

    # Get the parameters of this ellipse
    x_center = shape.coord_list[0]
    y_center = shape.coord_list[1]
    x_radius = shape.coord_list[2]
    y_radius = shape.coord_list[3]

    # Create a box to estimate the background
    x_min = int(round(x_center - x_radius))
    x_max = int(round(x_center + x_radius))
    y_min = int(round(y_center - y_radius))
    y_max = int(round(y_center + y_radius))

    # Return the extents
    return x_min, x_max, y_min, y_max

# *****************************************************************

def get_enclosing_boxes(region):

    """
    This function ...
    :param region:
    :return:
    """

    boxes = []

    # This is a hack to use mpl to determine the outer bounds of the regions
    # (but it's a legit hack - pyregion needs a major internal refactor before
    # we can approach this any other way)
    mpl_objs = region.get_mpl_patches_texts(origin=0)[0]

    # Loop over all objects
    for obj in mpl_objs:

        # Find the minimal enclosing box containing the shape
        extent = obj.get_extents()
        x_min, y_min = extent.min
        x_max, y_max = extent.max

        # Add the extent of this box
        boxes.append((x_min, x_max, y_min, y_max))

    return boxes

# *****************************************************************

def create_annulus(self, region, outer_factor, inner_factor=1.0):

    """
    This function ...
    :param region:
    :param outer_factor:
    :param inner_factor:
    :return:
    """

    # Create a new region
    region_annulus = pyregion.ShapeList([])

    # ...
    for shape in region:

        # Create new shapes by deep-copying the original shape
        # Creating new shapes from scratch: pyregion.parser_helper.Shape(None, None)
        inner_shape = copy.deepcopy(shape)
        outer_shape = copy.deepcopy(shape)

        # Add a '-' symbol to the name of the inner region
        inner_shape.name = '-' + shape.name

        # Set the size of the inner shape
        inner_shape.coord_list[2] *= inner_factor
        inner_shape.coord_list[3] *= inner_factor

        # Do special things to make this an excluded region (We're not supposed to do this as well)
        inner_shape.exclude = True

        # Set the size of the outer shape
        outer_shape.coord_list[2] *= outer_factor
        outer_shape.coord_list[3] *= outer_factor

        region_annulus.append(outer_shape)
        region_annulus.append(inner_shape)

    # Return the new region
    return region_annulus

# *****************************************************************

def expand(region, factor):

    """
    This function ...
    :param region:
    :param factor:
    :return:
    """

    # Create a new region
    region_expanded = pyregion.ShapeList([])

    # ...
    for shape in region:

        # Create a new shape
        expanded_shape = copy.deepcopy(shape)

        # Set the size of the new shape
        expanded_shape.coord_list[2] *= factor
        expanded_shape.coord_list[3] *= factor

        # Add the new shape to the new region
        region_expanded.append(expanded_shape)

    # Return the new region
    return region_expanded

# *****************************************************************

def ellipses_from_coordinates(coordinates):

    """
    This function creates a region consisting of ellipses, based on a list of coordinates
    :param coordinates: the list of coordinates for the different ellipses
    :return: the region of ellipses
    """

    # Initialize the region string
    region_string = "# Region file format: DS9 version 3.0\n"
    region_string += "global color=green\n"
    region_string += "image\n"

    # Loop over the objects in the coordinates list, adding a line for each one
    for object in coordinates:

        line = "ellipse(" + str(object.x_mean.value) + "," + str(object.y_mean.value) + "," + str(object.x_stddev.value) + "," + str(object.y_stddev.value) + ",0.0)\n"
        region_string += line

    # Parse the region string into a region object
    region = pyregion.parse(region_string)

    # Return the region
    return region

# *****************************************************************

def create_mask(region, header, x_size, y_size):

    """
    This function ...
    :param region:
    :param header:
    :param x_size:
    :param y_size:
    :return:
    """

    # Create a mask and return it
    return region.get_mask(header=header, shape=(y_size,x_size))

# *****************************************************************