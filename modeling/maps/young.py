#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.youngstars Contains the YoungStellarMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from .component import MapsComponent
from ...core.tools import filesystem as fs
from ...core.plot.distribution import DistributionPlotter
from ...magic.core.image import Image
from ...magic.maps.youngstars.young import YoungStellarMapsMaker

# -----------------------------------------------------------------

class YoungStellarMapMaker(MapsComponent):

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
        super(YoungStellarMapMaker, self).__init__(config, interactive)

        # -- Attributes --

        # The input FUV and FUV error maps
        self.fuv = None
        self.fuv_errors = None

        # The NORMALIZED (to unity) DISK IMAGE
        self.disk = None

        # The maps of the corrected FUV emission
        self.corrected_fuv_maps = dict()

        # The distributions of corrected FUV pixel values
        self.corrected_fuv_distributions = dict()

        # The map of young stars
        self.map = None

        # The path to the maps/young/fuv directory
        self.maps_young_fuv_path = None

        # Region of area taken for calculating distribution of pixel values
        self.distribution_region = None

        # The image of significance masks
        self.significance = Image()

        # The cutoff mask
        self.cutoff_mask = None

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

        # Load the necessary maps
        self.load_maps()

        # 3. Calculate the significance masks
        #self.calculate_significance()

        # 3. Make the map of young stars
        self.make_maps()

        # 5. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(YoungStellarMapMaker, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def load_frames(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the necessary data ...")

        # Load the GALEX FUV image and error map
        self.load_fuv()

        # Load the disk image and normalize to unity
        #self.load_disk()

    # -----------------------------------------------------------------

    def load_fuv(self):

        """
        This function ...
        :return:
        """

        # Get FUV frame and error map
        self.fuv = self.dataset.get_frame("GALEX FUV") # in original MJy/sr units
        self.fuv_errors = self.dataset.get_errormap("GALEX FUV") # in original MJy/sr units

    # -----------------------------------------------------------------

    def load_maps(self):

        """
        THis function ...
        :return:
        """

        # Inform the user
        log.info("Loading maps ...")

        # Load FUV attenuation map
        self.load_fuv_attenuation_map()

        # Load old stellar map
        self.load_old_stellar_map()

    # -----------------------------------------------------------------

    def load_fuv_attenuation_map(self):

        """
        THis function ...
        :return:
        """

        # Inform the user
        log.info("Loading the map of the FUV attenuation ...")

    # -----------------------------------------------------------------

    def load_old_stellar_map(self):

        """
        This function ...
        """

        # Inform the user
        log.info("Loading the map of old stars ...")

    # -----------------------------------------------------------------

    def calculate_significance(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the significance masks ...")

        # Get the significance mask
        if self.config.fuv_significance > 0: self.significance.add_mask(self.dataset.get_significance_mask("GALEX FUV", self.config.fuv_significance), "GALEX_FUV")

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Making the maps ...")

        maker = YoungStellarMapsMaker()

        maker.run()

        self.maps = maker.maps

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the maps
        self.write_maps()

        # Write the corrected FUV maps
        self.write_fuv_maps()

        # Write distribution region
        self.write_distribution_region()

        # Write histograms of corrected 24 micron pixels
        self.write_24mu_histograms()

        # Write the final young stellar map
        #self.write_map()

        # Write the significance mask
        self.write_significance_masks()

        # Write the cutoff mask
        self.write_cutoff_mask()

    # -----------------------------------------------------------------

    def write_fuv_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the corrected FUV maps ...")

        # Loop over the corrected FUV maps
        for factor in self.corrected_fuv_maps:

            # Determine the path
            path = fs.join(self.maps_young_fuv_path, str(factor) + ".fits")

            # Write
            self.corrected_fuv_maps[factor].saveto(path)

    # -----------------------------------------------------------------

    def write_distribution_region(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the distribution region ...")

        path = fs.join(self.maps_young_fuv_path, "histogram.reg")
        self.distribution_region.saveto(path)

    # -----------------------------------------------------------------

    def write_24mu_histograms(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the histograms of the corrected 24 micron pixels in the specified region ...")

        # Create a distribution plotter
        plotter = DistributionPlotter()

        # Loop over the distributions
        for factor in self.corrected_fuv_distributions:

            # Determine path
            path = fs.join(self.maps_young_fuv_path, str(factor) + " histogram.pdf")

            # Plot the distribution as a histogram
            plotter.add_distribution(self.corrected_fuv_distributions[factor], "Correction factor of " + str(factor))
            plotter.run(path)

            # Clear the distribution plotter
            plotter.clear()

    # -----------------------------------------------------------------

    def write_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the map of young stars ...")

        # Write
        self.map.saveto(self.young_stellar_map_path)

    # -----------------------------------------------------------------

    def write_significance_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the significance masks ...")

        # Write
        self.significance.saveto(self.young_stellar_significance_path)

    # -----------------------------------------------------------------

    def write_cutoff_mask(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the cutoff mask ...")

        # Write
        self.cutoff_mask.saveto(self.young_stellar_cutoff_path)

# -----------------------------------------------------------------
