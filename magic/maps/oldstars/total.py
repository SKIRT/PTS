#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.stars.old Contains the OldStellarMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....core.tools.logging import log
from ....magic.core.image import Image
from ....core.basics.configurable import Configurable

# -----------------------------------------------------------------

class TotalOldStellarMapMaker(Configurable):

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
        super(TotalOldStellarMapMaker, self).__init__(config, interactive)

        # -- Attributes --

        # The IRAC I1 frame in Jy
        self.i1_jy = None

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

        # 2. Load the necessary frames
        self.load_frames()

        # 3. Calculate the significance masks
        self.calculate_significance()

        # 4. Make the map of old stars
        self.make_map()

        # 5. Normalize the map
        self.normalize_map()

        # Make the cutoff mask
        self.make_cutoff_mask()

        # 6. Cut-off the map
        self.cutoff_map()

        # 7. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(TotalOldStellarMapMaker, self).setup(**kwargs)

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
        #conversion_factor = 1.0
        #conversion_factor *= 1e6

        # Convert the 3.6 micron image from Jy / sr to Jy / pixel
        #pixelscale = frame.average_pixelscale
        #pixel_factor = (1.0 / pixelscale ** 2).to("pix2/sr").value
        #conversion_factor /= pixel_factor

        # DO THE CONVERSION
        #frame *= conversion_factor
        #frame.unit = "Jy"

        # Set the frame
        self.i1_jy = frame.convert_to("Jy")

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

        self.map = self.i1_jy

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

    def make_cutoff_mask(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the cutoff mask ...")

        # Combine the significance masks
        high_significance = self.significance.intersect_masks()

        # Fill holes
        if self.config.remove_holes: high_significance.fill_holes()

        # Set
        self.cutoff_mask = high_significance.inverse()

    # -----------------------------------------------------------------

    def cutoff_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Cutting-off the map at low significance of the data ...")

        # Set zero outside of significant pixels
        self.map[self.cutoff_mask] = 0.0

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

        # Write the significance mask
        self.write_significance_masks()

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

    def write_significance_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the significance masks ...")

        # Write
        self.significance.saveto(self.old_stellar_significance_path)

# -----------------------------------------------------------------
