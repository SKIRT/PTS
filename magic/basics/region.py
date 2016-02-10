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
from astropy import coordinates as coord

# Import the relevant AstroMagic classes and modules
from .vector import Position, Extent
from .geometry import Ellipse, SkyEllipse
from .mask import Mask

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
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        raise ValueError("This function is not ready yet")

        # Create a new sky region
        region = cls()

        # Open the region file with pyregion and check if its in image coordinates
        _region = pyregion.open(path)
        def check_sky_coord(_region):
            if [s for s in _region if s.coord_format != "fk5"]: return False
            else: return True
        if not check_sky_coord(_region): raise ValueError("Region is not in sky coordinates")

        for shape in _region:

            x_center = shape.coord_list[0]
            y_center = shape.coord_list[1]
            x_radius = shape.coord_list[2]
            if shape.name == "ellipse":
                y_radius = shape.coord_list[3]
                try: angle = Angle(shape.coord_list[4], u.Unit("deg"))
                except: angle = Angle(0.0, u.Unit("deg"))
            elif shape.name == "circle":
                y_radius = shape.coord_list[2]
                angle = Angle(0.0, u.Unit("deg"))
            else: raise ValueError("Shapes other than ellipses or circles are not supported yet")

            center = coord.SkyCoord(ra=ra_center, dec=dec_center, unit=(u.Unit("deg"), u.Unit("deg")), frame='fk5')
            radius = Extent(x_radius, y_radius)

            # Only works for ellipses for now
            ellipse = Ellipse(center, radius, angle)

            # Add the ellipse to the region
            region.append(ellipse)

        # Return the new region
        return region

    # -----------------------------------------------------------------

    def append(self, shape):

        """
        This function ...
        :return:
        """

        # Check whether the shape is in sky coordinates
        if not isinstance(shape, SkyEllipse): raise ValueError("Shape must be of type SkyEllipse (for now)")

        # Otherwise, add the shape
        super(SkyRegion, self).append(shape)

    # -----------------------------------------------------------------

    def __mul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        new_region = SkyRegion()

        for ellipse in self: new_region.append(ellipse * value)

        # Return the new region
        return new_region

    # -----------------------------------------------------------------

    def __div__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        new_region = SkyRegion()

        for ellipse in self: new_region.append(ellipse / value)

        # Return the new region
        return new_region

    # -----------------------------------------------------------------

    def to_image_coordinates(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        # Initialize a new list contain the ellipses in image coordinates
        new_region = newRegion()

        # Fill the new list
        for i in range(len(self)): new_region.append(self[i].to_ellipse(wcs))

        # Return the list of ellipses in image coordinates
        return new_region

    # -----------------------------------------------------------------

    def to_mask(self, wcs):

        """
        This function ...
        :return:
        """

        x_size = wcs.naxis1
        y_size = wcs.naxis2

        mask = Mask(np.zeros((y_size, x_size)))

        for shape in self:

            ellipse = shape.to_ellipse(wcs)
            mask += Mask.from_ellipse(x_size, y_size, ellipse)

        # Return the mask
        return mask

# -----------------------------------------------------------------

class newRegion(list):

    """
    This class ...
    """

    def __init__(self):

        """
        This function ...
        :return:
        """

        # Call the constructor of the base class
        super(newRegion, self).__init__()

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Create a new region
        region = cls()

        # Open the region file with pyregion and check if its in image coordinates
        _region = pyregion.open(path)
        if not _region.check_imagecoord(): raise ValueError("Region is not in image coordinates")

        for shape in _region:

            # Skip shapes that are not ellipses or circles
            if shape.name != "ellipse" and shape.name != "circle": continue

            x_center = shape.coord_list[0]
            y_center = shape.coord_list[1]
            x_radius = shape.coord_list[2]
            if shape.name == "ellipse":
                y_radius = shape.coord_list[3]
                try: angle = Angle(shape.coord_list[4], u.Unit("deg"))
                except: angle = Angle(0.0, u.Unit("deg"))
            elif shape.name == "circle":
                y_radius = shape.coord_list[2]
                angle = Angle(0.0, u.Unit("deg"))
            else: raise ValueError("Shapes other than ellipses or circles are not supported yet")

            center = Position(x_center, y_center)
            radius = Extent(x_radius, y_radius)

            # Only works for ellipses for now
            ellipse = Ellipse(center, radius, angle)

            # Add the ellipse to the region
            region.append(ellipse)

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
        if not isinstance(shape, Ellipse): raise ValueError("Shape must be of type Ellipse (for now)")

        # Otherwise, add the shape
        super(newRegion, self).append(shape)

    # -----------------------------------------------------------------

    def __mul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        new_region = newRegion()

        for ellipse in self: new_region.append(ellipse * value)

        # Return the new region
        return new_region

    # -----------------------------------------------------------------

    def __div__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        new_region = newRegion()

        for ellipse in self: new_region.append(ellipse / value)

        # Return the new region
        return new_region

    # -----------------------------------------------------------------

    def to_sky_coordinates(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        # Initialize a new list to contain the ellipses in sky coordinates
        new_region = SkyRegion()

        # Fill the new list
        for i in range(len(self)): new_region.append(SkyEllipse.from_ellipse(self[i], wcs))

        # Return the list of ellipses in sky coordinates
        return new_region

    # -----------------------------------------------------------------

    def to_mask(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        x_size = wcs.naxis1
        y_size = wcs.naxis2

        mask = Mask(np.zeros((y_size, x_size)))

        for shape in self:

            mask += Mask.from_ellipse(x_size, y_size, shape)

        # Return the mask
        return mask

    # -----------------------------------------------------------------

    def save(self, path):

        """
        This function ...
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
            #suffix = " # "
            #color_suffix = "color = white"
            #text_suffix = "text = {" + text + "}"
            #suffix += color_suffix + " " + text_suffix
            print("image;ellipse({},{},{},{},{})".format(center.x+1, center.y+1, major, minor, angle), file=f)

        # Close the file
        f.close()

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
