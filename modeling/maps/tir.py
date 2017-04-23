#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.tir Contains the TIRMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import MapsComponent
from ...core.tools.logging import log
from ...magic.maps.tir.single import SingleBandTIRMapMaker
from ...magic.maps.tir.multi import MultiBandTIRMapMaker

# -----------------------------------------------------------------

class TIRMapMaker(MapsComponent):

    """
    This class...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(TIRMapMaker, self).__init__(config, interactive)

        # The TIR maps
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

        # 2. Make maps based on a single band
        self.make_maps_single()

        # 3. Make maps based on multiple bands
        self.make_maps_multi()

        # 4. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(TIRMapMaker, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def make_maps_single(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Making maps based on a single band ...")

        # Create
        maker = SingleBandTIRMapMaker()

        # Run
        maker.run(distance=self.galaxy_distance)

        # Set the maps
        self.maps["single"] = maker.maps

    # -----------------------------------------------------------------

    def make_maps_multi(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Making maps based on multiple bands ...")

        # Create
        maker = MultiBandTIRMapMaker()

        # Run
        maker.run(distance=self.galaxy_distance)

        # Set the maps
        self.maps["multi"] = maker.maps

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
