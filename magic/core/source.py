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
from .box import Box
from ..tools import plotting

# *****************************************************************

class Source(object):

    """
    This class...
    """

    def __init__(self, center, radius=None, x_radius=None, y_radius=None, angle=None, cutout=None,
                 background=None, background_mask=None, background_fit=None, source_mask=None, x_peak=None, y_peak=None):

        """
        The constructor ...
        :return:
        """

        # Set the necessary attributes
        self.center = center

        # Set attribute for circle
        self.radius = radius

        # Set attributes for ellipses
        self.x_radius = x_radius
        self.y_radius = y_radius
        self.angle = angle

        # Set other attributes
        self.cutout = cutout
        self.background = background
        self.background_mask = background_mask
        self.background_fit = background_fit
        self.source_mask = source_mask
        self.x_peak = x_peak
        self.y_peak = y_peak

        self.removed = None

    # *****************************************************************

    def find_peak(self):

        """
        This function ...
        :return:
        """



    # *****************************************************************

    def plot(self):

        """
        This function ...
        :return:
        """

        cutout_background = self.background_fit[self.cutout.y_min-self.background.y_min:self.cutout.y_max-self.background.y_min,
                                         self.cutout.x_min-self.background.x_min:self.cutout.x_max-self.background.x_min]

        # Calculate the relative position of the center for the cutout and background boxes
        rel_center = self.cutout.rel_position(self.center)
        rel_center_background = self.background.rel_position(self.center)

        # Do the plotting
        plotting.plot_source(rel_center.x, rel_center.y, rel_center_background.x, rel_center_background.y, self.background,
                             self.background_mask, self.background_fit, self.cutout, cutout_background, self.source_mask, self.removed)

# *****************************************************************