#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.dust.cortese Contains the AttenuationDustMapsMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy import constants

# Import the relevant PTS classes and modules
from ....magic.core.frame import Frame
from ....core.tools import introspection, tables
from ....core.tools import filesystem as fs
from ....core.tools.logging import log
from ....magic.tools.colours import make_colour_map
from ....core.units.parsing import parse_unit as u
from ....core.basics.configurable import Configurable

# -----------------------------------------------------------------

speed_of_light = constants.c
solar_luminosity = 3.846e26 * u("W")

# -----------------------------------------------------------------

#ssfr_colours = ["FUV-H", "FUV-i", "FUV-r", "FUV-g", "FUV-B"]
ssfr_colours = ["FUV-H", "FUV-i", "FUV-r", "FUV-g"]

colour_combinations = {"FUV-H": ("GALEX FUV", "2MASS H"),
                       "FUV-i": ("GALEX FUV", "SDSS i"),
                       "FUV-r": ("GALEX FUV", "SDSS r"),
                       "FUV-g": ("GALEX FUV", "SDSS g")}

# -----------------------------------------------------------------

def make_map():

    """
    This function ...
    :return: 
    """

    maker = AttenuationDustMapsMaker()

    maker.run()

# -----------------------------------------------------------------

class AttenuationDustMapsMaker(Configurable):

    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(AttenuationDustMapsMaker, self).__init__(*args, **kwargs)

        # -- Attributes --

        # Frames and error maps
        self.frames = dict()
        self.errors = dict()

        # The maps
        self.maps = dict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # Make the maps
        self.make_maps()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs: 
        :return: 
        """

        # Call the setup function of the base class
        super(AttenuationDustMapsMaker, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Making the maps ...")

# -----------------------------------------------------------------
