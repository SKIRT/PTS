#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.basics.properties Contains the GalaxyProperties class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.composite import SimplePropertyComposite

# -----------------------------------------------------------------

class GalaxyProperties(SimplePropertyComposite):
    
    """
    This class ...
    """

    def __init__(self, **kwargs):

        """
        The constructor ...
        """

        self.name = kwargs.pop("name", None)
        self.ngc_id = kwargs.pop("ngc_id", None)
        self.center = kwargs.pop("center", None)
        self.major = kwargs.pop("major", None)
        self.major_arcsec = kwargs.pop("major_arcsec", None)
        self.ellipticity = kwargs.pop("ellipticity", None)
        self.position_angle = kwargs.pop("position_angle", None)
        self.distance = kwargs.pop("distance", None)
        self.distance_error = kwargs.pop("distance_error", None)
        self.inclination = kwargs.pop("inclination", None)

# -----------------------------------------------------------------
