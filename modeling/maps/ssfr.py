#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.ssfr Contains the SSFRMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from .component import MapsComponent
from ...magic.maps.ssfr.colours import ColoursSSFRMapsMaker

# -----------------------------------------------------------------

class SSFRMapMaker(MapsComponent):

    """
    This class...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(SSFRMapMaker, self).__init__(config, interactive)

        # THe maps
        self.maps = dict()

    # -----------------------------------------------------------------

    @property
    def maps_sub_path(self):

        """
        This function ...
        :return: 
        """

        return self.maps_ssfr_path

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Make SSFR maps based on colours
        self.make_ssfr_colours()

        # 3. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(SSFRMapMaker, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def make_ssfr_colours(self):

        """
        THis function ...
        :return: 
        """

        # Inform the user
        log.info("Making sSFR maps based on colours ...")

        # Create the map maker
        maker = ColoursSSFRMapsMaker()

        # Run the maker
        maker.run()

        # Get the maps
        self.maps = maker.maps

    # -----------------------------------------------------------------

    def write(self):

        """
        THis function ...
        :return: 
        """

        # Inform the user
        log.info("Writing ...")

        # Write the maps
        self.write_maps()

        # Write the origins
        self.write_origins()

# -----------------------------------------------------------------
