#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.stars.old Contains the OldStellarMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....core.tools.logging import log
from ..component import MapsComponent
from .disk import DiskOldStellarMapMaker
from .bulge import BulgeOldStellarMapMaker
from .total import TotalOldStellarMapMaker

# -----------------------------------------------------------------

class OldStellarMapMaker(MapsComponent):

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
        super(OldStellarMapMaker, self).__init__(config, interactive)

        # The maps
        self.maps =

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # Make disk map
        self.make_disk_map()

        # Make total map
        self.make_total_map()

        # Make bulge map
        self.make_bulge_map()

        # 7. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(OldStellarMapMaker, self).setup()

    # -----------------------------------------------------------------

    def make_disk_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the old stellar disk map ...")

        # Create the maker
        maker = DiskOldStellarMapMaker()

        # Run
        maker.run()

    # -----------------------------------------------------------------

    def make_total_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the total map of the old stars ...")

        # Create the maker
        maker = TotalOldStellarMapMaker()

        # Run
        maker.run()

    # -----------------------------------------------------------------

    def make_bulge_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the map of the old stellar bulge ...")

        # Create the maker
        maker = BulgeOldStellarMapMaker()

        # Run
        maker.run()

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

# -----------------------------------------------------------------
