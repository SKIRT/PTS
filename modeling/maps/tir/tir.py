#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.tir.tir Contains the TIRMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..component import MapsComponent
from ....core.tools.logging import log
from .single import SingleBandTIRMapMaker
from .multi import MultiBandTIRMapMaker

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

        # -- Attributes --

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
        maker.run()

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
