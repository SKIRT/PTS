#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.dust.attenuation Contains the AttenuationDustMapsMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....core.tools.logging import log
from ....core.basics.configurable import Configurable

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

        self.attenuation = None
        self.attenuation_origins = None

        # The maps
        self.maps = dict()

        # The origins
        self.origins = dict()

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

        # Get input
        self.attenuation = kwargs.pop("attenuation")

        # Get origins
        self.attenuation_origins = kwargs.pop("attenuation_origins", None)

        # Get already created maps
        self.maps = kwargs.pop("maps", dict())

    # -----------------------------------------------------------------

    @property
    def has_origins(self):

        """
        This function ...
        :return:
        """

        return self.attenuation_origins is not None

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Making the dust maps ...")

        # Loop over the colour maps
        for name in self.attenuation:

            # Set origin
            if self.has_origins:

                # Set the origins
                origins = self.attenuation_origins[name]
                self.origins[name] = origins

            # Check whether a dust map is already present
            if name in self.maps:
                log.warning("The " + name + " dust map is already created: not creating again")
                continue

            # Get the map
            dust = self.attenuation[name]

            # Normalized
            dust.normalize()

            # Set as dust map
            self.maps[name] = dust

# -----------------------------------------------------------------
