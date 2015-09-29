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
from .trackrecord import TrackRecord

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

        # Initialize a track record of sources
        self.track_record = None

    # *****************************************************************

    @property
    def has_source(self):

        """
        This function ...
        :return:
        """

        return self.source is not None

    # *****************************************************************

    @property
    def has_model(self):

        """
        This function ...
        :return:
        """

        return self.model is not None

    # *****************************************************************

    @property
    def fwhm(self):

        """
        This function ...
        :return:
        """

        pass

    # *****************************************************************

    def enable_track_record(self):

        """
        This function ...
        :return:
        """

        # Create a new track record
        self.track_record = TrackRecord()

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
        self.source = analysis.find_source(frame, center, radius, 0.0, config, self.track_record)

    # *****************************************************************

    def fit_model(self, config):

        """
        This function ...
        :param frame:
        :param config:
        :return:
        """

        # Fit model to the source, in a loop over different analytical forms for the model
        for level in range(len(config.model_names)):

            # Do the fitting
            source, model = analysis.fit_model_to_source(self.source, config, self.track_record, level=level)

            # If a model was found, set the attributes of the star object and exit the loop
            if model is not None:

                self.source = source
                self.model = model
                break

    # *****************************************************************

    def remove(self, frame):

        """
        This function removes the star from a given frame
        :param frame:
        :return:
        """

        pass

# *****************************************************************