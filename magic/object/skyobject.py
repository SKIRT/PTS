#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.sky.skyobject Contains the abstract SkyObject class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta
from abc import abstractmethod

# Import the relevant PTS classes and modules
from ..basics.stretch import PixelStretch
from ..basics.trackrecord import TrackRecord
from ..analysis import sources

# -----------------------------------------------------------------

class SkyObject(object):

    """
    This class...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, position):

        """
        The constructor ...
        """

        # Position
        self.position = position

        # Set the source attribute to None initially
        self.source = None

        # Set the contour (for saturation) to None initially
        self.contour = None

        # Initialize a track record of sources
        self.track_record = None

        # Set the special flag to False initially
        self.special = False

        # Set the ignore flag to False initially
        self.ignore = False

    # -----------------------------------------------------------------

    @property
    def has_source(self):

        """
        This function ...
        :return:
        """

        return self.source is not None

    # -----------------------------------------------------------------

    @property
    def has_contour(self):

        """
        This function ...
        :return:
        """

        return self.contour is not None

    # -----------------------------------------------------------------

    @property
    def has_track_record(self):

        """
        This function ...
        :return:
        """

        return self.track_record is not None

    # -----------------------------------------------------------------

    def enable_track_record(self):

        """
        This function ...
        :return:
        """

        # Create a new track record
        self.track_record = TrackRecord("Start")

    # -----------------------------------------------------------------

    @abstractmethod
    def ellipse(self, wcs, initial_radius):

        """
        This function ...
        :param wcs:
        :param initial_radius:
        :return:
        """

        return

    # -----------------------------------------------------------------

    @abstractmethod
    def ellipse_parameters(self, wcs, initial_radius):

        """
        This function ...
        :param wcs:
        :param initial_radius:
        :return:
        """

        return

    # -----------------------------------------------------------------

    def pixel_position(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        # Get the x and y coordinate of the object's position
        return self.position.to_pixel(wcs)

    # -----------------------------------------------------------------

    def find_source(self, frame, config):

        """
        This function ...
        :param frame:
        :param config:
        :return:
        """

        # Add a new stage to the track record
        if self.has_track_record: self.track_record.set_stage("detection")

        # Get the parameters of the circle
        radius = PixelStretch(config.initial_radius, config.initial_radius)
        ellipse = self.ellipse(frame.wcs, radius)

        # Find a source
        self.source = sources.find_source(frame, ellipse, config, self.track_record, special=self.special)

    # -----------------------------------------------------------------

    def find_contour(self, frame, config, saturation=False):

        """
        This function ...
        :param frame:
        :param config:
        :param saturation:
        :return:
        """

        # Determine which box and mask
        if saturation:

            box = self.saturation.cutout
            mask = self.saturation.mask
        else:

            box = self.source.cutout
            mask = self.source.mask

        # Get the aperture
        contour = sources.find_contour(box, mask, config.sigma_level)

        if contour is None: return

        # Calculate the difference (in number of pixels) between the aperture center and the position of the sky object
        difference = self.pixel_position(frame.wcs) - contour.center

        # Set the aperture if the offset is smaller than or equal to the specified maximum
        if difference.norm <= config.max_offset or config.max_offset is None:

            # Set the aperture
            self.contour = contour

# -----------------------------------------------------------------
