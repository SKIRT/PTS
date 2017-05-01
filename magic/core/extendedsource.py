#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.core.extendedsource Contains the ExtendedSource class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.coordinates import Angle

# Import the relevant PTS classes and modules
from .source import Source
from .detection import Detection
from ..basics.stretch import PixelStretch
from ..region.ellipse import PixelEllipseRegion

# -----------------------------------------------------------------

class ExtendedSource(Source):

    """
    This class...
    """

    def __init__(self, **kwargs):

        """
        The constructor ...
        :param position:
        """

        # Call the constructor of the base class
        super(ExtendedSource, self).__init__(**kwargs)

        # Set other attributes
        self.name = kwargs.pop("name", None)
        self.redshift = kwargs.pop("redshift", None)
        self.galaxy_type = kwargs.pop("galaxy_type", None)
        self.names = kwargs.pop("names", None)
        self.distance = kwargs.pop("distance", None)
        self.inclination = kwargs.pop("inclination", None)
        self.d25 = kwargs.pop("d25", None)
        self.major = kwargs.pop("major", None)
        self.minor = kwargs.pop("minor", None)
        self.pa = kwargs.pop("position_angle", None)
        self.principal = kwargs.pop("principal", False)
        self.companions = kwargs.pop("companions", [])
        self.parent = kwargs.pop("parent", None)

    # -----------------------------------------------------------------

    @property
    def companion(self):

        """
        This function ...
        :return:
        """

        return self.parent is not None

    # -----------------------------------------------------------------

    def pa_for_wcs(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        try: orientation = wcs.standard_orientation_angle
        except ValueError: orientation = wcs.orientation_angle

        # Add the orientation angle (w.r.t. standard E-W and S-N projection on the x and y axes) to the position angle
        # that is expressed in the standard way
        return self.pa + orientation

    # -----------------------------------------------------------------

    def contains(self, position):

        """
        This function ...
        :param position:
        :return:
        """

        return self.detection.contains(position)

    # -----------------------------------------------------------------

    @property
    def has_extent(self):

        """
        This function ...
        :return:
        """

        # Check whether the length of the major axis is defined
        return self.major is not None

    # -----------------------------------------------------------------

    def ellipse(self, wcs, default_radius):

        """
        This function ...
        :param wcs:
        :param default_radius:
        :return:
        """

        center, radius, angle = self.ellipse_parameters(wcs, default_radius)
        return PixelEllipseRegion(center, radius, angle)

    # -----------------------------------------------------------------

    @property
    def shape(self):

        """
        This property ...
        :return:
        """

        return self.detection.shape

    # -----------------------------------------------------------------

    def ellipse_parameters(self, wcs, default_radius):

        """
        This function ...
        :param wcs:
        :param default_radius:
        :return:
        """

        if self.pa is None: angle = Angle(0.0, "deg")
        else: angle = self.pa_for_wcs(wcs)

        if self.major is None:

            x_radius = default_radius
            y_radius = default_radius

        elif self.minor is None or angle == 0.0:

            x_radius = 0.5 * self.major.to("arcsec").value / wcs.pixelscale.x.to("arcsec").value
            y_radius = x_radius

        else:

            x_radius = 0.5 * self.major.to("arcsec").value / wcs.pixelscale.x.to("arcsec").value
            y_radius = 0.5 * self.minor.to("arcsec").value / wcs.pixelscale.y.to("arcsec").value

        pixel_position = self.pixel_position(wcs)

        # Return the parameters
        return pixel_position, PixelStretch(x=x_radius, y=y_radius), angle

    # -----------------------------------------------------------------

    def detection_from_shape(self, frame, shape, outer_factor):

        """
        This function ...
        :param frame:
        :param shape:
        :param outer_factor:
        :return:
        """

        # Create the source
        self.detection = Detection.from_shape(frame, shape, outer_factor)

    # -----------------------------------------------------------------

    def detection_from_parameters(self, frame, outer_factor, expansion_factor=1.0):

        """
        This function ...
        :param frame:
        :param outer_factor:
        :param expansion_factor:
        :return:
        """

        # Get the parameters describing the elliptical contour
        ellipse = self.ellipse(frame.wcs, None)

        #print(ellipse, dir(ellipse))
        #print(ellipse.center)
        #print(ellipse.radius)
        #print(ellipse.angle)
        #print(self.has_extent)
        #print(self.major)

        if ellipse.center.x < 0 or ellipse.center.y < 0:
            self.detection = None
            return

        # Create a source object
        ellipse *= expansion_factor
        self.detection = Detection.from_ellipse(frame, ellipse, outer_factor)

    # -----------------------------------------------------------------

    def remove(self, frame, mask, config):

        """
        This function ...
        :param frame:
        :param mask:
        :param config:
        :return:
        """

        # If a segment was found that can be identified with a source
        if self.has_detection or config.remove_if_undetected:

            # Estimate the background
            self.detection.estimate_background(config.interpolation_method, config.sigma_clip)

            # Replace the frame with the estimated background
            self.detection.background.replace(frame, where=self.detection.mask)

            # Update the mask
            mask[self.detection.cutout.y_slice, self.detection.cutout.x_slice] += self.detection.mask

# -----------------------------------------------------------------
