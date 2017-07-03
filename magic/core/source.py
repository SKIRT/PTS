#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.core.source Contains the Source class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta
from abc import abstractmethod

# Import the relevant PTS classes and modules
from ..basics.stretch import PixelStretch
from ..analysis import sources

# -----------------------------------------------------------------

class Source(object):

    """
    This class...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Set attributes
        self.index = kwargs.pop("index")
        self.position = kwargs.pop("position")
        self.special = kwargs.pop("special", False)
        self.ignore = kwargs.pop("ignore", False)
        
        # The detection
        self.detection = None

        # The contour
        self.contour = None

    # -----------------------------------------------------------------

    @property
    def has_detection(self):

        """
        This function ...
        :return:
        """

        return self.detection is not None

    # -----------------------------------------------------------------

    @property
    def has_contour(self):

        """
        This function ...
        :return:
        """

        return self.contour is not None

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

    def detect(self, frame, config):

        """
        This function ...
        :param frame:
        :param config:
        :return:
        """

        # Add a new stage to the track record
        #if self.has_track_record: self.track_record.set_stage("detection")
        track_record = None

        # Get the parameters of the circle
        radius = PixelStretch(config.initial_radius, config.initial_radius)
        ellipse = self.ellipse(frame.wcs, radius)

        # Find a source
        self.detection = sources.find_source(frame, ellipse, config, track_record, special=self.special)

    # -----------------------------------------------------------------

    def find_contour(self, frame, config):

        """
        This function ...
        :param frame:
        :param config:
        :return:
        """

        # Box and mask
        box = self.detection.cutout
        mask = self.detection.mask

        # Implementation
        self._find_contour_impl(frame.wcs, box, mask, config)

    # -----------------------------------------------------------------

    def _find_contour_impl(self, wcs, box, mask, config):

        """
        This function ...
        :param wcs:
        :param box:
        :param mask:
        :param config:
        :return:
        """

        # Get the aperture
        contour = sources.find_contour(box, mask, config.sigma_level)

        # No contour found
        if contour is None: return

        # Calculate the difference (in number of pixels) between the aperture center and the position of the sky object
        difference = self.pixel_position(wcs) - contour.center

        # Set the aperture if the offset is smaller than or equal to the specified maximum
        if difference.norm <= config.max_offset or config.max_offset is None:

            # Set the aperture
            self.contour = contour

# -----------------------------------------------------------------
