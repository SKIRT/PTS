#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Import Python 3 functionality
from __future__ import (absolute_import, division, print_function)

# Import standard modules
import numpy as np

# Import astronomical units
from astropy import units as u

# Import Astromagic modules
from ..tools import analysis
from .vector import Position, Extent
from ..tools import interpolation

# *****************************************************************

class Galaxy(object):

    """
    This class ...
    """

    def __init__(self, pgc_id=None, center=None, names=None, major=None, minor=None, pa=None):

        """
        The constructor ...
        :param ra:
        :param dec:
        :param name:
        :return:
        """

        # Set the attributes
        self.pgc_id = pgc_id
        self.center = center
        self.names = names
        self.major = major
        self.minor = minor
        self.pa = pa
        self.principal = False

        # Set the source attribute to None initially
        self.source = None

        # Set the model attribute to None initially
        self.model = None

    # *****************************************************************

    def ellipse_parameters(self, wcs, pixelscale, default_radius):

        """
        This function ...
        :param default_radius:
        :return:
        """

        # Get the center of the galaxy in pixel coordinates
        x_center, y_center = self.center.to_pixel(wcs, origin=0)

        if self.pa is None: angle = 0.0
        else: angle = self.pa.value

        if self.major is None:

            x_radius = default_radius
            y_radius = default_radius

        elif self.minor is None or self.pa == 0.0:

            print(pixelscale)

            x_radius = self.major.to("arcsec") / pixelscale
            y_radius = x_radius

        else:

            x_radius = self.major.to("arcsec") / pixelscale
            y_radius = self.minor.to("arcsec") / pixelscale

        # Return the parameters
        return Position(x=x_center, y=y_center), Extent(x=x_radius, y=y_radius), angle

    # *****************************************************************

    def find_source(self, frame, config):

        """
        This function ...
        :return:
        """

        # Get the parameters of the ellipse
        center, radius, angle = self.ellipse_parameters(frame.wcs, frame.pixelscale, config.initial_radius)

        # Find a source
        self.source = analysis.find_source(frame, center, radius, angle, config)

    # *****************************************************************

    def fit_model(self, config):

        """
        This function ...
        :param frame:
        :param config:
        :return:
        """

        pass

    # *****************************************************************

    def remove(self, frame, config):

        """
        This function ...
        :param frame:
        :param config:
        :return:
        """

        subtraction_method = config.subtraction_method
        subtract_if_undetected = config.subtract_if_undetected

        if self.source is None and subtract_if_undetected:

            # Get the parameters of the ellipse
            x_center, y_center, x_radius, y_radius, angle = self.ellipse_parameters(frame.wcs, frame.pixelscale, config.initial_radius)

            #self.source = Source(x_center, y_center, x_radius, y_radius, angle, cutout=box, background=background_box, background_mask=background_mask,
            #                     background_fit=polynomial, source_mask=box_mask)

        elif self.source is not None:

            # Calculate the interpolated background
            interpolated_cutout = interpolation.in_paint(self.source.cutout, self.source.source_mask)

            self.source.removed = interpolated_cutout

            frame[self.source.cutout.y_min:self.source.cutout.y_max,self.source.cutout.x_min:self.source.cutout.x_max] = interpolated_cutout

    # *****************************************************************

    def plot(self):
        
        pass
        
# *****************************************************************