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
from astropy.coordinates import Angle

# Import Astromagic modules
from .source import Source
from ..tools import analysis
from .vector import Position
from ..tools import fitting
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
    def has_track_record(self):

        """
        This function ...
        :return:
        """

        return self.track_record is not None

    # *****************************************************************

    @property
    def fwhm(self):

        """
        This function ...
        :return:
        """

        # Return the fwhm value of the model
        return fitting.fwhm(self.model)

    # *****************************************************************

    @property
    def flux(self):

        """
        This function ...
        :return:
        """

        # Return the flux of the source
        return self.source.flux

    # *****************************************************************

    def enable_track_record(self):

        """
        This function ...
        :return:
        """

        # Create a new track record
        self.track_record = TrackRecord("Start")

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

        # Add a new stage to the track record
        if self.has_track_record: self.track_record.set_stage("detection")

        # Get the parameters of the circle
        center, radius = self.circle_parameters(frame.wcs, frame.pixelscale, config.initial_radius)

        # Find a source
        self.source = analysis.find_source(frame, center, radius, Angle(0.0, u.deg), config, self.track_record, special=False)

    # *****************************************************************

    def fit_model(self, config):

        """
        This function ...
        :param frame:
        :param config:
        :return:
        """

        # Add a new stage to the track record
        if self.has_track_record: self.track_record.set_stage("fitting")

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

    def remove(self, frame, config, default_fwhm):

        """
        This function removes the star from a given frame
        :param frame:
        :return:
        """

        # Add a new stage to the track record
        if self.has_track_record: self.track_record.set_stage("removal")

        # Convert FWHM to sigma
        default_sigma = default_fwhm / 2.35

        radius = fitting.sigma(self.model) * config.sigma_level if self.model is not None else default_sigma * config.sigma_level

        # Determine the center position of the source (center of model if present, otherwise position of the star)
        if self.source is not None:

            # If the star has been modeled succesfully, use the center position of the model
            # Otherwise, use the source's peak
            if self.model is not None: center = fitting.center(self.model)
            else: center = self.source.peak

        else:

            # Calculate the pixel coordinate of the star's position
            position_x, position_y = self.position.to_pixel(frame.wcs, origin=0)
            center = Position(x=position_x, y=position_y)

        # Create a source
        source = Source(frame, center, radius, Angle(0.0, u.deg), config.outer_factor)

        # Estimate the background
        source.estimate_background(config.method, config.sigma_clip)

        # Add the source to the track record
        if self.has_track_record: self.track_record.append(source)

        # Replace the frame with the estimated background
        source.estimated_background.replace(frame, where=source.background_mask)

        # Use the new source
        self.source = source

    # *****************************************************************

    def remove_saturation(self, frame, config, default_fwhm):

        """
        This function ...
        """

        # Convert FWHM to sigma
        default_sigma = default_fwhm / 2.35

        # Determine the radius for the saturation detection
        model = self.model
        radius = fitting.sigma(model) * config.sigmas if model is not None else default_sigma * config.sigmas

        # Add a new stage to the track record
        if self.has_track_record: self.track_record.set_stage("saturation")

        # Look for a center segment corresponding to a 'saturation' source
        source = analysis.find_source_segmentation(frame, self.source.center, radius, Angle(0.0, u.deg), config, track_record=self.track_record)

        # If a 'saturation' source was found
        if source is not None:

            # Replace the source by a source that covers the saturation
            self.source = source

            # Estimate the background
            self.source.estimate_background(config.remove_method, config.sigma_clip)

            # Replace the frame with the estimated background
            self.source.estimated_background_cutout.replace(frame, where=self.source.mask)

# *****************************************************************