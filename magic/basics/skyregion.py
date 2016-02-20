#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.basics.skyregion Contains the SkyRegion class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
import pyregion
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord

# Import the relevant AstroMagic classes and modules
from .vector import Extent, Position
from .skygeometry import SkyLine, SkyCircle, SkyEllipse, SkyRectangle, SkyPolygon

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

            # Get the coordinate format of this shape
            coord_format = shape.coord_format

            # The shape is a point -> SkyCoord
            if shape.name == "point":

                # Get RA and declination
                ra = shape.coord_list[0]
                dec = shape.coord_list[1]

                # Create sky coordinate and add it to the region
                coordinate = SkyCoord(ra=ra, dec=dec, unit="deg", frame=coord_format)
                region.append(coordinate)

            # The shape is a line or vector -> SkyLine
            elif shape.name == "line" or shape.name == "vector":

                # Get the RA and declination of the two points
                ra_1 = shape.coord_list[0]
                dec_1 = shape.coord_list[1]
                ra_2 = shape.coord_list[2]
                dec_2 = shape.coord_list[3]

                # Create the sky coordinates
                coordinate_1 = SkyCoord(ra=ra_1, dec=dec_1, unit="deg", frame=coord_format)
                coordinate_2 = SkyCoord(ra=ra_2, dec=dec_2, unit="deg", frame=coord_format)

                # Create the SkyLine object and add it to the region
                line = SkyLine(coordinate_1, coordinate_2)
                region.append(line)

            # The shape is a circle -> SkyCircle
            elif shape.name == "circle":

                # Get the RA and declination of the center and the radius
                ra_center = shape.coord_list[0]
                dec_center = shape.coord_list[1]
                radius = shape.coord_list[2] * u.Unit("deg")

                # Create a sky cooridnate for the center
                center = SkyCoord(ra=ra_center, dec=dec_center, unit="deg", frame=coord_format)

                # Create a SkyCircle object and add it to the region
                circle = SkyCircle(center, radius)
                region.append(circle)

            # The shape is an ellipse -> SkyEllipse
            elif shape.name == "ellipse":

                # Get the RA and declination of the center
                ra_center = shape.coord_list[0]
                dec_center = shape.coord_list[1]
                center = SkyCoord(ra=ra_center, dec=dec_center, unit="deg", frame=coord_format)

                # Get the radius
                x_radius = shape.coord_list[2] * u.Unit("deg")
                y_radius = shape.coord_list[3] * u.Unit("deg")
                radius = Extent(x_radius, y_radius)

                # Get the angle
                angle = Angle(shape.coord_list[4], "deg")

                # Create a SkyEllipse object and add it to the region
                ellipse = SkyEllipse(center, radius, angle)
                region.append(ellipse)

            # The shape is a rectangle -> SkyRectangle
            elif shape.name == "box":

                # Get the RA and declination of the center
                ra_center = shape.coord_list[0]
                dec_center = shape.coord_list[1]
                center = SkyCoord(ra=ra_center, dec=dec_center, unit="deg", frame=coord_format)

                # Get the width and height
                width = shape.coord_list[2] * u.Unit("deg")
                height = shape.coord_list[3] * u.Unit("deg")

                # Create radius
                radius = Extent(0.5 * width, 0.5*height)

                # Create a SkyRectangle and add it to the region
                rectangle = SkyRectangle(center, radius)
                region.append(rectangle)

            # The shape is a polygon -> SkyPolygon
            elif shape.name == "polygon":

                number_of_points = 0.5 * len(shape.coord_list)
                assert int(number_of_points) == number_of_points
                number_of_points = int(number_of_points)

                # Create a new SkyPolygon
                polygon = SkyPolygon()

                # Get the RA and declination of the different points
                for i in range(number_of_points):

                    # Create a new SkyCoord object
                    ra = shape.coord_list[0]
                    dec = shape.coord_list[1]
                    coordinate = SkyCoord(ra=ra, dec=dec, unit="deg", frame=coord_format)

                    # Add the coordinate to the polygon
                    polygon.add_point(coordinate)

                # Add the polygon to the region
                region.append(polygon)

            # Unrecognized shape
            else: raise ValueError("Unrecognized shape (should be point, line, vector, circle, ellipse, box or polygon")

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
        if not shape.__class__.__name__.startswith("Sky"): raise ValueError("Shape must be SkyCoord, SkyLine, SkyCircle, SkyEllipse or SkyRectangle")

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

    def to_image_coordinates(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        # Avoid circular import at module level
        from .region import Region

        # Initialize a new list contain the ellipses in image coordinates
        new_region = Region()

        # Fill the new list
        for shape in self:

            if shape.__class__.__name__ == "SkyCoord":
                x, y = shape.to_pixel(wcs, origin=0, mode='wcs')
                new_region.append(Position(x, y))
            elif shape.__class__.__name__ == "SkyLine": new_region.append(shape.to_line(wcs))
            elif shape.__class__.__name__ == "SkyCircle": new_region.append(shape.to_circle(wcs))
            elif shape.__class__.__name__ == "SkyEllipse": new_region.append(shape.to_ellipse(wcs))
            elif shape.__class__.__name__ == "SkyRectangle": new_region.append(shape.to_rectangle(wcs))
            elif shape.__class__.__name__ == "SkyPolygon": new_region.append(shape.to_polygon(wcs))
            else: raise ValueError("Uncrecognized shape")

        # Return the list of ellipses in image coordinates
        return new_region

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
            ra_deg = ellipse.center.ra.to("deg").value
            dec_deg = ellipse.center.dec.to("deg").value
            major = ellipse.major.to("arcsec").value
            minor = ellipse.minor.to("arcsec").value
            angle = ellipse.angle.degree

            line = "fk5;ellipse(%s,%s,%.2f\",%.2f\",%s)" % (ra_deg, dec_deg, major, minor, angle)

            # Write to region file
            print(line, file=f)

        # Close the file
        f.close()

# -----------------------------------------------------------------
