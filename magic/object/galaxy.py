#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.sky.galaxy Contains the Galaxy class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .skyobject import SkyObject

# -----------------------------------------------------------------

class Galaxy(SkyObject):

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
        super(Galaxy, self).__init__(**kwargs)

        # Set properties
        self.name = kwargs.pop("name", None)
        self.redshift = kwargs.pop("redshift", None)
        self.type = kwargs.pop("galaxy_type", None)
        self.names = kwargs.pop("names", None)
        self.distance = kwargs.pop("distance", None)
        self.inclination = kwargs.pop("inclination", None)
        self.d25 = kwargs.pop("d25", None)
        self.major = kwargs.pop("major", None)
        self.minor = kwargs.pop("minor", None)
        self.pa = kwargs.pop("position_angle", None)
        self.principal = kwargs.pop("principal", False)
        self.companions = kwargs.pop("companions", [])
        self.parent = kwargs.pop("parent", None)

    # -----------------------------------------------------------------

    @property
    def companion(self):

        """
        This function ...
        :return:
        """

        return self.parent is not None

    # -----------------------------------------------------------------

    def pa_for_wcs(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        try: orientation = wcs.standard_orientation_angle
        except ValueError: orientation = wcs.orientation_angle

        # Add the orientation angle (w.r.t. standard E-W and S-N projection on the x and y axes) to the position angle
        # that is expressed in the standard way
        return self.pa + orientation

# -----------------------------------------------------------------
