#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.maps.ssfr Contains the SSFRMapsAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import MapsAnalysisComponent
from ....core.basics.log import log
from ....magic.maps.ssfr.colours import ColoursSSFRMapsMaker

# -----------------------------------------------------------------

class SSFRMapsAnalyser(MapsAnalysisComponent):

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
        super(SSFRMapsAnalyser, self).__init__(*args, **kwargs)

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
        super(SSFRMapsAnalyser, self).setup(**kwargs)

        # Load the analysis run
        self.load_run()

    # -----------------------------------------------------------------

    def make_ssfr_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making sSFR maps ...")

        # Set the method name
        method_name = "ssfr"

        # Get current
        if self.config.remake: current = dict()
        else: current = self.current_maps

        # Get origins and methods
        colours_origins = self.get_colours_origins(flatten=True)
        colours_methods = self.get_colours_methods(flatten=True)

        # Create the map maker
        maker = ColoursSSFRMapsMaker()

        # Run the maker
        maker.run(colours=self.colours, colours_origins=colours_origins, colours_methods=colours_methods, maps=current, method_name=method_name)

        # Get the maps
        self.maps = maker.maps

        # Get the origins
        self.origins = maker.origins

        # Set the methods
        self.methods = maker.methods

# -----------------------------------------------------------------
