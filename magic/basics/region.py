#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import astronomical modules
import pyregion
import astropy.coordinates as coord

# Import AstroMagic modules
from .vector import Extent

# -----------------------------------------------------------------

class Region(pyregion.ShapeList):

    """
    This class ...
    """

    # -----------------------------------------------------------------

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

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, filepath, wcs=None):

        """
        This function ...
        :param filepath:
        :return:
        """

        region = pyregion.open(filepath)

        # Convert to image coordinates if a WCS is given
        if wcs: region = region.as_imagecoord(wcs.to_header())

        # Change the class name of the ShapeList instance to our own Region class
        region.__class__ = cls

        # Set the attributes that are added in the Region subclass
        region.selected = False

        # Return the Region object
        return region

    # -----------------------------------------------------------------

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

    # -----------------------------------------------------------------

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

    # -----------------------------------------------------------------

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

    # -----------------------------------------------------------------

    def select(self):

        """
        This function ...
        :return:
        """

        self.selected = True

    # -----------------------------------------------------------------

    def deselect(self):

        """
        This function ...
        :return:
        """

        self.selected = False

    # -----------------------------------------------------------------

    def expand(self, factor):

        """
        This function ...
        """

        pass

# -----------------------------------------------------------------
