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

        # Call the constructor of the base class
        super(GalaxyProperties, self).__init__()

        # Define properties
        self.add_property("name", "string", "name of the galaxy", None)
        self.add_property("ngc_name", "string", "NGC name of the galaxy", None)
        self.add_property("hyperleda_name", "string", "HYPERLEDA name of the galaxy", None)
        self.add_property("galaxy_type", "string", "galaxy type", None)
        self.add_property("center", "skycoordinate", "center coordinate", None)
        self.add_property("major", "quantity", "major axis length", None)
        self.add_property("major_arcsec", "quantity", "major axis length (arcseconds)", None)
        self.add_property("ellipticity", "real", "ellipticity", None)
        self.add_property("position_angle", "angle", "position angle", None)
        self.add_property("distance", "quantity", "distance", None)
        self.add_property("distance_error", "quantity", "distance error", None)
        self.add_property("inclination", "angle", "inclination", None)
        self.add_property("redshift", "real", "redshift", None)
        self.add_property("common_name", "string", "common name", None)

        # Set values
        self.set_properties(kwargs)

# -----------------------------------------------------------------
