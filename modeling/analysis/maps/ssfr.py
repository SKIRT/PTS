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
from ....magic.maps.ssfr.colours import ColoursSSFRMapsMaker, ssfr_colours

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

        # The colour maps
        self.colours = dict()

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        Thisf ucntion ...
        :param kwargs:
        :return:
        """

        # 2. Load the colour maps
        self.load_colours()

        # 3. Make the maps
        self.make_maps()

        # Analyse
        self.analyse_maps()

        # 4. Write
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(SSFRMapsAnalyser, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def load_colours(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the colour maps ...")

        # Loop over the possible colours for tracing sSFR
        for colour in ssfr_colours:

            # Get colour map and the name
            colour_map, colour_name = self.get_colour_map_and_name(colour)

            # Check if found
            if colour_map is None:
                log.warning("Could not find a '" + colour + "' colour map")
                continue

            # Add the colour map
            self.colours[colour_name] = colour_map

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making sSFR maps ...")

        # Set the method name
        method_name = "colours"

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

    def analyse_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Analysing the maps ...")

    # -----------------------------------------------------------------

    @property
    def maps_sub_path(self):

        """
        This function ...
        :return:
        """

        return self.ssfr_path

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the colour maps
        self.write_maps()

        # Write origins
        self.write_origins()

        # Write the methods
        self.write_methods()

# -----------------------------------------------------------------
