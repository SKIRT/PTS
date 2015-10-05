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

    @property
    def has_source(self):

        """
        This function ...
        :return:
        """

        return self.source is not None

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

        if self.source is None and config.subtract_if_undetected:

            # Get the parameters of the ellipse
            center, radius, angle = self.ellipse_parameters(frame.wcs, frame.pixelscale, config.default_radius)

            pass

        # If a segment was found that can be identified with a source
        if self.has_source:

            # Estimate the background
            self.source.estimate_background(config.remove_method, config.sigma_clip)

            # Replace the frame with the estimated background
            self.source.estimated_background_cutout.replace(frame, where=self.source.mask)

    # *****************************************************************

    def plot(self):

        """
        This function ...
        :return:
        """

        pass
        
# *****************************************************************