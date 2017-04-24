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

# Import standard modules
import numpy as np

# Import astronomical modules
from astropy import constants

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from .component import MapsComponent
from ...core.tools import filesystem as fs
from ...core.basics.distribution import Distribution
from ...core.plot.distribution import DistributionPlotter
from ...magic.core.image import Image
from ...core.units.parsing import parse_unit as u
from ...magic.maps.ionizingstars.ionizing import IonizingStellarMapsMaker

# -----------------------------------------------------------------

speed_of_light = constants.c
solar_luminosity = 3.846e26 * u("W")

# -----------------------------------------------------------------

class IonizingStellarMapMaker(MapsComponent):

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
        super(IonizingStellarMapMaker, self).__init__(config, interactive)

        # -- Attributes --

        # The Halpha map IN SOLAR UNITS
        self.halpha = None

        # The maps of hot dust
        self.hots = None

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

        # Calculate the significance masks
        #self.calculate_significance()

        # 3. Make the map
        self.make_maps()

        # ...
        #self.create_distribution_region()
        #self.make_distributions()

        # 4. Normalize the map
        #self.normalize_map()

        # Make the cutoff mask
        #self.make_cutoff_mask()

        # 5. Cut-off map
        #self.cutoff_map()

        # 5. Writing
        self.write()

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

    def load_frames(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the necessary data ...")

        # Load the MIPS 24 micron image and convert to solar units
        self.load_hot()

        # Load the H alpha image and convert to solar units
        self.load_halpha()

    # -----------------------------------------------------------------

    def calculate_significance(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the significance masks ...")

        # Get the significance mask
        if self.config.mips24_significance > 0: self.significance.add_mask(self.dataset.get_significance_mask("MIPS 24mu", self.config.mips24_significance), "MIPS_24mu")
        if self.config.halpha_significance > 0: self.significance.add_mask(self.get_halpha_significance_mask(self.config.halpha_significance), "Halpha")

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

    # -----------------------------------------------------------------

    def load_halpha(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the H-alpha image and converting to solar units ...")

        # Get the H-alpha image
        self.halpha = self.masked_halpha_frame

        # Convert from erg/s to Lsun
        self.halpha.convert_to("Lsun")

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return: 
        """

        # Create
        maker = IonizingStellarMapsMaker()

        # Run
        maker.run(halpha=self.halpha, hots=self.hots)

        # Set the maps
        self.maps = maker.maps

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        self.write_maps()

# -----------------------------------------------------------------
