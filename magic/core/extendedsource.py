#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.core.extendedsource Contains the ExtendedSource class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .source import Source

# -----------------------------------------------------------------

class ExtendedSource(Source):

    """
    This class...
    """

    def __init__(self, **kwargs):

        """
        The constructor ...
        :param position:
        """

        # Call the constructor of the base class
        super(ExtendedSource, self).__init__(**kwargs)

        # Set other attributes
        self.name = kwargs.pop("name", None)
        self.redshift = kwargs.pop("redshift", None)
        self.galaxy_type = kwargs.pop("galaxy_type", None)
        self.names = kwargs.pop("names", None)
        self.distance = kwargs.pop("distance", None)
        self.inclination = kwargs.pop("inclination", None)
        self.d25 = kwargs.pop("d25", None)
        self.major = kwargs.pop("major", None)
        self.minor = kwargs.pop("minor", None)
        self.position_angle = kwargs.pop("position_angle", None)
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
