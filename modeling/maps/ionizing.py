#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.stars.ionizing Contains the IonizingStellarMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from .component import MapsComponent
from ...magic.maps.ionizingstars.ionizing import IonizingStellarMapsMaker
from ...core.filter.filter import parse_filter

# -----------------------------------------------------------------

methods = None

# -----------------------------------------------------------------

class IonizingStellarMapMaker(MapsComponent):

    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(IonizingStellarMapMaker, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The Halpha map, if present
        self.halpha = None

        # The maps of hot dust
        self.hots = None

        # Origins
        self.hots_origins = None

        # Methods
        self.hots_methods = None

    # -----------------------------------------------------------------

    @property
    def maps_sub_path(self):

        """
        This function ...
        :return: 
        """

        return self.maps_ionizing_path

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Load the necessary frames
        self.load_frames()

        # 3. Make the map
        self.make_maps()

        # 5. Writing
        self.write()

        # 6. Plotting
        if self.config.plot: self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(IonizingStellarMapMaker, self).setup(**kwargs)

    # -----------------------------------------------------------------

    @property
    def has_halpha(self):

        """
        This function ...
        :return:
        """

        return self.has_frame_for_filter("Halpha")

    # -----------------------------------------------------------------

    def load_frames(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the necessary data ...")

        # Load the MIPS 24 micron image (and convert to solar units == > NO?)
        self.load_hot()

        # Load the H alpha image (and convert to solar units == > NO?)
        if self.has_halpha: self.load_halpha()

    # -----------------------------------------------------------------

    def load_hot(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Loading the maps of hot dust ...")

        # Get
        self.hots = self.get_hot_dust_maps()

        # Get
        self.hots_origins = self.get_hot_dust_origins()

        # Methods
        self.hots_methods = self.get_hot_dust_methods()

    # -----------------------------------------------------------------

    def load_halpha(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the H-alpha image and converting to solar units ...")

        # Get the H-alpha image
        self.halpha = self.get_frame_for_filter(parse_filter("Halpha"))

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Making the map ...")

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

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write
        self.write_maps()

        # Write origins
        self.write_origins()

        # Write the methods
        self.write_methods()

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot the maps
        self.plot_maps(scale="linear")

        # Plot the contours
        self.plot_contours(filled=True)

        # Plot the radial profiles
        self.plot_profiles()

# -----------------------------------------------------------------
