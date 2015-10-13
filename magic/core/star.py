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
from ..tools import statistics
from .skyobject import SkyObject
from .source import Source
from ..tools import analysis
from .vector import Position
from ..tools import fitting

# *****************************************************************

class Star(SkyObject):

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
        self.position_error = position_error
        self.ra_error = ra_error
        self.dec_error = dec_error
        self.k_mag = k_mag
        self.b_mag = b_mag
        self.v_mag = v_mag
        self.r_mag = r_mag
        self.i_mag = i_mag

        # Set the model attribute to None initially
        self.model = None

        # Set the has_saturation flag to False initially
        self.has_saturation = False

        # Call the constructor of the base class
        super(Star, self).__init__(position)

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

    def ellipse_parameters(self, wcs, pixelscale, default_radius):

        """
        This function ...
        :param wcs:
        :param pixelscale:
        :param initial_radius:
        :return:
        """

        # Return the parameters
        return self.pixel_position(wcs), default_radius, Angle(0.0, u.deg)

    # *****************************************************************

    def fit_model(self, config, source=None):

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

            if source is None:

                # Do the fitting
                source, model = analysis.fit_model_to_source(self.source, config, self.track_record, level=level)

            else:

                source, model = analysis.fit_model_to_source(source, config, self.track_record, level=level)

            # If a model was found, set the attributes of the star object and exit the loop
            if model is not None:

                self.source = source
                self.model = model
                break

    # *****************************************************************

    def remove(self, frame, config, default_fwhm, sigma_clip):

        """
        This function removes the star from a given frame
        :param frame:
        :return:
        """

        # If a segment was found that can be identified with a source
        if self.has_source or config.remove_if_undetected:

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
                center = self.pixel_position(frame.wcs)

            # Create a source
            source = Source(frame, center, radius, Angle(0.0, u.deg), config.outer_factor)

            # Estimate the background
            source.estimate_background(config.method, sigma_clip)

            # Add the source to the track record
            if self.has_track_record: self.track_record.append(source)

            # Replace the frame with the estimated background
            source.background.replace(frame, where=source.mask)

    # *****************************************************************

    def remove_saturation(self, frame, config, default_fwhm, sigma_clip):

        """
        This function ...
        """

        # Convert FWHM to sigma
        default_sigma = default_fwhm * statistics.fwhm_to_sigma

        # Determine the radius for the saturation detection
        model = self.model
        radius = fitting.sigma(model) * config.sigmas if model is not None else default_sigma * config.sigmas

        # Add a new stage to the track record
        if self.has_track_record: self.track_record.set_stage("saturation")

        # Look for a center segment corresponding to a 'saturation' source
        source = analysis.find_source_segmentation(frame, self.pixel_position(frame.wcs), radius, Angle(0.0, u.deg), config, track_record=self.track_record, special=self.special)

        # If a 'saturation' source was found
        if source is not None:

            # Replace the source by a source that covers the saturation
            self.source = source

            # Estimate the background
            self.source.estimate_background(config.remove_method, sigma_clip)

            # Replace the frame with the estimated background
            self.source.background.replace(frame, where=self.source.mask)

            # Indicate that a saturation source was found
            return True

        # Otherwise, return False
        else: return False

# *****************************************************************