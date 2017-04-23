#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.maps.oldstars.bulge Contains the BulgeOldStellarMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....core.tools.logging import log
from ....magic.core.image import Image
from ....core.basics.configurable import Configurable

# -----------------------------------------------------------------

class BulgeOldStellarMapMaker(Configurable):

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
        super(BulgeOldStellarMapMaker, self).__init__(config, interactive)

        # -- Attributes --

        # The map of the old stars
        self.map = None

        # The image of significance masks
        self.significance = Image()

        # The cutoff mask
        self.cutoff_mask = None

        # THe maps
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

        # 4. Make the map of old stars
        self.make_map()

        # 5. Normalize the map
        self.normalize_map()

        # Make the cutoff mask
        #self.make_cutoff_mask()

        # 6. Cut-off the map
        #self.cutoff_map()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(BulgeOldStellarMapMaker, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def load_frames(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the necessary data ...")

        # Load the IRAC I1 frame and convert to Jansky

        frame = self.dataset.get_frame("IRAC I1")

        # Convert the 3.6 micron image from MJy/sr to Jy/sr
        conversion_factor = 1.0
        conversion_factor *= 1e6

        # Convert the 3.6 micron image from Jy / sr to Jy / pixel
        pixelscale = frame.average_pixelscale
        pixel_factor = (1.0 / pixelscale ** 2).to("pix2/sr").value
        conversion_factor /= pixel_factor

        # DO THE CONVERSION
        frame *= conversion_factor
        frame.unit = "Jy"

        # Set the frame
        self.i1_jy = frame

    # -----------------------------------------------------------------

    def calculate_significance(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the significance masks ...")

        # Get the significance mask
        if self.config.i1_significance > 0: self.significance.add_mask(self.dataset.get_significance_mask("IRAC I1", self.config.i1_significance), "IRAC_I1")

    # -----------------------------------------------------------------

    def make_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the map of old stars ...")

        # Create copy
        self.map = self.masked_bulge_frame

    # -----------------------------------------------------------------

    def normalize_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Normalizing the map of old stars ...")

        # Normalize the old stellar map
        self.map.normalize()
        self.map.unit = None

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the map of old stars
        self.write_map()

    # -----------------------------------------------------------------

    def write_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the map of old stars ...")

        # Write
        #self.map.saveto(self.old_stellar_map_path)

# -----------------------------------------------------------------
