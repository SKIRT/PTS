#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.sky.star Contains the abstract Star class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .skyobject import SkyObject

# -----------------------------------------------------------------

class Star(SkyObject):

    """
    This class ...
    """

    def __init__(self, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(Star, self).__init__(**kwargs)

        # Set the attributes
        self.catalog = kwargs.pop("catalog", None)
        self.id = kwargs.pop("id", None)
        self.ra_error = kwargs.pop("ra_error", None)
        self.dec_error = kwargs.pop("dec_error", None)

        # The FWHM table
        self.fwhms = kwargs.pop("fwhms", None)

        #self.magnitudes = magnitudes
        #self.magnitude_errors = magnitude_errors

        #self.on_galaxy = on_galaxy

# -----------------------------------------------------------------
