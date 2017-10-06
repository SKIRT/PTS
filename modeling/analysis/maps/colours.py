#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.maps.colours Contains the ColourMapsAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import MapsAnalysisComponent
from ....core.basics.log import log

# -----------------------------------------------------------------

class ColourMapsAnalyser(MapsAnalysisComponent):

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
        super(ColourMapsAnalyser, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        Thisf ucntion ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ColourMapsAnalyser, self).setup(**kwargs)

        # Load the analysis run
        self.load_run()

    # -----------------------------------------------------------------

    def make_colour_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making colour maps ...")

        # Current maps
        if self.config.remake: current = dict()
        else: current = self.current_maps

        # Create the map maker
        maker = ColourMapsMaker()

        # Run the map maker
        maker.run(colours=self.available_colours, frames=self.frames, maps=current)

        # Set the maps
        self.maps = maker.maps

        # Set the origins
        self.origins = maker.origins

        # Set the methods
        self.methods = maker.methods

# -----------------------------------------------------------------
