#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import MapsAnalysisComponent
from ....core.basics.log import log

# -----------------------------------------------------------------

class IonizingMapsAnalyser(MapsAnalysisComponent):

    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(IonizingMapsAnalyser, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        Thisf ucntion ...
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
        super(IonizingMapsAnalyser, self).setup(**kwargs)

        # Load the analysis run
        self.load_run()

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Making ionizing stellar maps ...")

        # Get the current maps
        if self.config.remake: current = dict()
        else: current = self.current_maps

        # Create
        maker = IonizingStellarMapsMaker()

        # Run
        maker.run(halpha=self.halpha, hots=self.hots, hots_origins=self.hots_origins, hots_methods=self.hots_methods, maps=current)

        # Set the maps
        self.maps = maker.maps

        # Set the origins
        self.origins = maker.origins

        # Set the methods
        self.methods = maker.methods

# -----------------------------------------------------------------
