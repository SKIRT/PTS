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
from .vector import Extent, Position
from .region import Region
from .geometry import Ellipse
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

        # Create a new sky region
        region = cls()

        # Open the region file with pyregion and check if its in image coordinates
        _region = pyregion.open(path)
        def check_sky_coord(_region):
            if [s for s in _region if s.coord_format != "fk5"]: return False
            else: return True
        if not check_sky_coord(_region): raise ValueError("Region is not in sky coordinates")

        for shape in _region:

            if shape.name == "point":

                ra = shape.coord_list[0]
                dec = shape.coord_list[1]

                coordinate = SkyCoord(ra=x, dec=y, unit=())

                region.append()

            if shape.name == "line" or shape.name == "vector":

                ra_1 = shape.coord_list[0]
                dec_1 = shape.coord_list[1]
                ra_2 = shape.coord_list[0]
                dec_2 = shape.coord_list[1]

            if shape.name == "circle":

                x_center = shape.coord_list[0]
                y_center = shape.coord_list[1]
                radius = shape.coord_list[2]

                angle = Angle(0.0, u.Unit("deg"))

            if shape.name == "ellipse":

                x_center = shape.coord_list[0]
                y_center = shape.coord_list[1]
                x_radius = shape.coord_list[2]
                y_radius = shape.coord_list[3]

                try: angle = Angle(shape.coord_list[4], u.Unit("deg"))
                except: angle = Angle(0.0, u.Unit("deg"))

            else: raise ValueError("Unrecognized shape")

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
            else: raise ValueError("Uncrecognized shape")

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
