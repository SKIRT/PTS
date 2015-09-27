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
from .source import Source
from ..tools import analysis
from .vector import Position

# *****************************************************************

class Star(object):

    """
    This class ...
    """

    def __init__(self, ucac_id=None, position=None, position_error=None, ra_error=None, dec_error=None, k_mag=None,
                 b_mag=None, v_mag=None, r_mag=None, i_mag=None):

        """
        The constructor ...
        :return:
        """

        # Set the attributes
        self.ucac_id = ucac_id
        self.position = position
        self.position_error = position_error
        self.ra_error = ra_error
        self.dec_error = dec_error
        self.k_mag = k_mag
        self.b_mag = b_mag
        self.v_mag = v_mag
        self.r_mag = r_mag
        self.i_mag = i_mag

        # Set the source attribute to None initially
        self.source = None

        # Set the model attribute to None initially
        self.model = None

    # *****************************************************************

    def circle_parameters(self, wcs, pixelscale, default_radius):

        """
        This function ...
        :param wcs:
        :param pixelscale:
        :param initial_radius:
        :return:
        """

        # Get the center of the galaxy in pixel coordinates
        x_center, y_center = self.position.to_pixel(wcs, origin=0)

        # Return the parameters
        return Position(x=x_center, y=y_center), default_radius

    # *****************************************************************

    def find_source(self, frame, config):

        """
        This function ...
        :return:
        """

        # Get the parameters of the circle
        center, radius = self.circle_parameters(frame.wcs, frame.pixelscale, config.initial_radius)

        # Find a source
        self.source = analysis.find_source(frame, center, radius, 0.0, config)

    # *****************************************************************

    def fit_model(self, config):

        """
        This function ...
        :param frame:
        :param config:
        :return:
        """

        # Fit a model
        self.model = analysis.fit_model_to_source(self.source, config)

# *****************************************************************