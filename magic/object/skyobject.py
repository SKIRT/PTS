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

# -----------------------------------------------------------------

class SkyObject(object):

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

        # Set properties
        self.index = kwargs.pop("index")
        self.position = kwargs.pop("position")
        self.sed = kwargs.pop("sed")

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
