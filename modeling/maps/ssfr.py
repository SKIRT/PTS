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
from ...magic.maps.ssfr.colours import ColoursSSFRMapMaker
from ...core.tools import filesystem as fs

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

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(SSFRMapMaker, self).setup()

    # -----------------------------------------------------------------

    def make_ssfr_colours(self):

        """
        THis function ...
        :return: 
        """

        # Inform the user
        log.info("Making sSFR maps based on colours ...")

        # Create the map maker
        maker = ColoursSSFRMapMaker()

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

    # -----------------------------------------------------------------

    def write_maps(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Writing the maps ...")

        # Loop over the maps
        for name in self.maps:

            # Determine path
            path = fs.join(self.maps_ssfr_path, name + ".fits")

            # Save
            self.maps[name].saveto(path)

# -----------------------------------------------------------------
