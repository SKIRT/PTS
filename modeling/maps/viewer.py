#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.viewer Contains the MapsViewer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from .component import MapsComponent
from ...magic.view.multi import MultiImageViewer

# -----------------------------------------------------------------

ncolumns = 3
image_width = 300
image_height = 300
background_color = "white"
key_color = "#4180d3"

# -----------------------------------------------------------------

class MapsViewer(MapsComponent):

    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(MapsViewer, self).__init__(*args, **kwargs)

        # Maps sub path
        self._maps_sub_path = None

        # The viewer
        self.viewer = None

    # -----------------------------------------------------------------

    @property
    def maps_sub_path(self):

        """
        This function ...
        :return: 
        """

        return self._maps_sub_path

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # Load the maps
        self.load_maps()

        # Generate the page
        self.generate_page()

        # Show
        self.show()

    # -----------------------------------------------------------------

    @property
    def has_colour_maps(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.has_colour_maps

    # -----------------------------------------------------------------

    @property
    def has_ssfr_maps(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.has_ssfr_maps

    # -----------------------------------------------------------------

    @property
    def has_tir_maps(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.has_tir_maps

    # -----------------------------------------------------------------

    @property
    def has_attenuation_maps(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.has_attenuation_maps

    # -----------------------------------------------------------------

    @property
    def has_old_maps(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.has_old_maps

    # -----------------------------------------------------------------

    @property
    def has_young_maps(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.has_young_maps

    # -----------------------------------------------------------------

    @property
    def has_ionizing_maps(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.has_ionizing_maps

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(MapsViewer, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Reset
        self._maps_sub_path = None

    # -----------------------------------------------------------------

    def load_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the maps ...")

    # -----------------------------------------------------------------

    def generate_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Viewing ...")

        viewer = MultiImageViewer()

        viewer.run()

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

# -----------------------------------------------------------------
