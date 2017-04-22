#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.ssfr.ssfr Contains the SSFRMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....core.tools.logging import log
from ..component import MapsComponent
from .colours import ColoursSSFRMapMaker

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
        if self.config.write: self.write()

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

        # Set the path
        maker.config.path = self.config.path

        # Run the maker
        maker.run()

    # -----------------------------------------------------------------

    def write(self):

        """
        THis function ...
        :return: 
        """

        # Inform the user
        log.info("Writing ...")

# -----------------------------------------------------------------
