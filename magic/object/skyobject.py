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
