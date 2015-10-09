#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Import standard modules
import copy
import numpy as np

# Import astronomical modules
import pyregion
import astropy.coordinates as coord

# Import Astromagic modules
from ..core.vector import Extent

# *****************************************************************

class Region(pyregion.ShapeList):

    """
    This class ...
    """

    # *****************************************************************

    def __init__(self):

        """
        The constructor ...
        :param region:
        :return:
        """

        # Call the ShapeList constructor
        super(Region, self).__init__()

        # Set as unselected initially
        self.selected = False

    # *****************************************************************

    @classmethod
    def ellipse(cls, center, radius, angle):

        """
        This function ...
        :param center:
        :param radius:
        :param angle:
        :return:
        """

        x_radius = radius.x if isinstance(radius, Extent) else radius
        y_radius = radius.y if isinstance(radius, Extent) else radius

        # Create a string identifying this ellipse
        region_string = "# Region file format: DS9 version 3.0\n"
        region_string += "global color=green\n"
        region_string += "image\n"
        region_string += "ellipse(" + str(center.x) + "," + str(center.y) + "," + str(x_radius) + "," + str(y_radius) + "," + str(angle.degree) + ")\n"

        # TODO: FIX THIS!! DOES RETURN A PYREGION OBJECT BUT WE WANT A ASTROMAGIC REGION OBJECT RETURNED !!

        # Create a region and return it
        return pyregion.parse(region_string)

    # *****************************************************************

    @classmethod
    def ellipses(cls, centers, radii, angles):

        """
        This function ...
        :param centers:
        :param radii:
        :param angles:
        :return:
        """

        # TODO: accept input in fk5 and create a region in sky coordinates

        # Create a string identifying this ellipse
        region_string = "# Region file format: DS9 version 3.0\n"
        region_string += "image\n"

        # Loop over the parameter sets
        for center, radius, angle in zip(centers, radii, angles):

            x_radius = radius.x if isinstance(radius, Extent) else radius
            y_radius = radius.y if isinstance(radius, Extent) else radius
            region_string += "ellipse(" + str(center.x) + "," + str(center.y) + "," + str(x_radius) + "," + str(y_radius) + "," + str(angle.degree) + ")\n"

        # Create a region and return it
        return pyregion.parse(region_string)

    # *****************************************************************

    @classmethod
    def circles(cls, centers, radii, colors=None):

        """
        This function ...
        :param centers:
        :param radii:
        :return:
        """

        # Check type of input
        if isinstance(centers[0], coord.SkyCoord): fk5 = True
        else: fk5 = False

        # Create a string identifying this ellipse
        region_string = "# Region file format: DS9 version 3.0\n"
        region_string += "global color=green\n"

        if fk5: region_string += "fk5\n"
        else: region_string += "image\n"

        if colors is None: colors = [colors] * len(centers)

        # Loop over the parameter sets
        for center, radius, color in zip(centers, radii, colors):

            if color is None: suffix = " # color = " + color
            else: suffix = ""

            if fk5: line = "fk5;circle({0},{1},{2:.2f}\")".format(center.ra.value, center.dec.value, radius.value) + suffix + "\n"
            else: line = "circle(" + str(center.x) + "," + str(center.y) + "," + str(radius) + ")" + suffix + "\n"

            #print(line)

            region_string += line

        # Create a region and return it
        return pyregion.parse(region_string)

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

def ellipse_parameters(shape):

    """
    This function ...
    :param shape:
    :return:
    """

    x_center = shape.coord_list[0]
    y_center = shape.coord_list[1]
    x_radius = shape.coord_list[2]

    if shape.name == "ellipse": y_radius = shape.coord_list[3]
    elif shape.name == "circle": y_radius = shape.coord_list[2]
    else: raise ValueError("Shape must be either a circle or an ellipse")
    
    return x_center, y_center, x_radius, y_radius

# *****************************************************************

def get_enclosing_box(shape):

    """
    This function ...
    :param shape:
    :return:
    """

    # TODO: make it work for shapes other than ellipses!

    # Get the parameters of this ellipse (or circle)
    x_center, y_center, x_radius, y_radius = ellipse_parameters(shape)

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

def create_annulus(region, outer_factor, inner_factor=1.0):

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

    # Loop over all shapes in the original region
    for shape in region:

        # Create a new shape
        expanded_shape = copy.deepcopy(shape)

        # Set the size of the new shape
        expanded_shape.coord_list[2] *= factor
        if shape.name == "ellipse": expanded_shape.coord_list[3] *= factor

        # Add the new shape to the new region
        region_expanded.append(expanded_shape)

    # Return the new region
    return region_expanded

# *****************************************************************

def ellipses(ra_list, dec_list, height_list, width_list, angle_list):

    # Initialize the region string
    region_string = "# Region file format: DS9 version 3.0\n"
    region_string += "global color=green\n"
    region_string += "fk5\n"

    for ra, dec, height, width, angle in zip(ra_list, dec_list, height_list, width_list, angle_list):

        line = "fk5;ellipse(%s,%s,%.2f\",%.2f\",%s)\n" % (ra, dec, height, width, angle)
        region_string += line

    region = pyregion.parse(region_string)

    # Return the region
    return region

# *****************************************************************

def circles(ra_list, dec_list, radius_list):

    # Initialize the region string
    region_string = "# Region file format: DS9 version 3.0\n"
    region_string += "global color=green\n"
    region_string += "fk5\n"

    for ra, dec, radius in zip(ra_list, dec_list, radius_list):

        line = "fk5;circle(%s,%s,%.2f\")\n" % (ra, dec, radius)
        region_string += line

    region = pyregion.parse(region_string)

    # Return the region
    return region

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

        if type(object).__name__ == "Gaussian2D": line = "ellipse(" + str(object.x_mean.value) + "," + str(object.y_mean.value) + "," + str(object.x_stddev.value) + "," + str(object.y_stddev.value) + ",0.0)\n"
        elif type(object).__name__ == "AiryDisk2D":
            # see https://en.wikipedia.org/wiki/Airy_disk#Approximation_using_a_Gaussian_profile and
            # http://astropy.readthedocs.org/en/latest/api/astropy.modeling.functional_models.AiryDisk2D.html#astropy.modeling.functional_models.AiryDisk2D
            sigma = 0.42 * object.radius.value * 0.81989397882
            line = "ellipse(" + str(object.x_0.value) + "," + str(object.y_0.value) + "," + str(sigma) + "," + str(sigma) + ",0.0)\n"
        else: raise ValueError("Models other than Gaussian2D and AiryDisk2D are not yet supported")

        region_string += line

    # Parse the region string into a region object
    region = pyregion.parse(region_string)

    # Return the region
    return region

# *****************************************************************

def one_ellipse(parameters):

    """
    This function ...
    :param parameters:
    :return:
    """

    # Create a string identifying this ellipse
    region_string = "# Region file format: DS9 version 3.0\n"
    region_string += "global color=green\n"
    region_string += "image\n"
    region_string += "ellipse(" + str(parameters[0]) + "," + str(parameters[1]) + "," + str(parameters[2]) + "," + str(parameters[3]) + "," + str(parameters[4]) + ")\n"

    # Create a region and return it
    return pyregion.parse(region_string)

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

def parse(region_string):

    """
    This function is a simple wrapper around the pyregion.parse function, to contain
    :param region_string:
    :return:
    """

    # Parse the region string and create a region
    return pyregion.parse(region_string)

# *****************************************************************

def scale(shape, factor):

    """
    This function ...
    :param shape:
    :param factor:
    :return:
    """

    new_shape = copy.deepcopy(shape)
    new_shape.coord_list[2] *= factor

    if new_shape.name == "ellipse": new_shape.coord_list[3] *= factor

    return new_shape

# *****************************************************************

def scale_circle(shape, factor):

    """
    This function ...
    :param shape:
    :param factor:
    :return:
    """

    new_shape = copy.deepcopy(shape)
    new_shape.coord_list[2] *= factor

    return new_shape

# *****************************************************************

def subtract(region_a, region_b, center_offset_tolerance, header):

    """
    This function ...
    :param region_a:
    :param region_b:
    :return:
    """

    # TODO: fix this function: do not only use the first shape of region b!!

    new_region = region_a.as_imagecoord(header)
    region_b = region_b.as_imagecoord(header)

    x_b = region_b[0].coord_list[0]
    y_b = region_b[0].coord_list[1]

    for i in range(len(new_region)):

        x_center = new_region[i].coord_list[0]
        y_center = new_region[i].coord_list[1]

        diff_x = x_center - x_b
        diff_y = y_center - y_b

        distance = np.sqrt(diff_x**2 + diff_y**2)

        if distance < center_offset_tolerance:

            del new_region[i]
            break

    # Return the subtracted region
    return new_region

# *****************************************************************

def in_box(shape, dimensions, require_completely=False):

    """
    This function ...
    :param shape:
    :param dimensions:
    :param require_completely:
    :return:
    """

    # TODO: use require_completely (require the shape to completely fall within the box)

    x_min, x_max, y_min, y_max = get_enclosing_box(shape)

    if x_min >= dimensions[1] or y_min >= dimensions[0]: return False
    if x_max < 0 or y_max < 0: return False

    return True

# *****************************************************************

def mean_radius(region):

    # Initialize an empty list to contain the different sigma values
    sigmas = []

    # Loop over all shapes in the region
    for shape in region:

        sigma_x = shape.coord_list[2]
        sigma_y = shape.coord_list[3]

        # Add the sigma, averaged over the x and y directions, to the list of sigmas
        sigmas.append(0.5*(sigma_x + sigma_y))

    return np.mean(sigmas)

# *****************************************************************

def max_radius(region):

    # Initialize an empty list to contain the different sigma values
    sigmas = []

    # Loop over all shapes in the region
    for shape in region:

        sigma_x = shape.coord_list[2]
        sigma_y = shape.coord_list[3]

        # Add the sigma, averaged over the x and y directions, to the list of sigmas
        sigmas.append(0.5*(sigma_x + sigma_y))

    return max(sigmas)

# *****************************************************************
