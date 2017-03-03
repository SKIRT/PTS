#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.stars.young Contains the YoungStellarMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....core.tools.logging import log
from ..component import MapsComponent
from ....core.tools import filesystem as fs
from ....core.basics.distribution import Distribution
from ....core.plot.distribution import DistributionPlotter
from ....magic.region.composite import PixelCompositeRegion
from ....magic.region.list import PixelRegionList
from ....magic.core.image import Image

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

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the necessary frames
        self.load_frames()

        # 3. Calculate the significance masks
        self.calculate_significance()

        # 3. Make the map of young stars
        self.make_map()

        # ...
        self.create_distribution_region()
        self.make_distributions()

        # 4. Normalize the map
        #self.normalize_map()

        # Create the cutoff mask
        self.make_cutoff_mask()

        # 5. Cut-off map
        #self.cutoff_map()

        # 5. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(YoungStellarMapMaker, self).setup()

        # Create
        self.maps_young_fuv_path = fs.create_directory_in(self.maps_young_path, "fuv")

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
        self.load_disk()

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

    def load_disk(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the disk image ...")

        # Get disk frame
        self.disk = self.masked_disk_frame

        # Normalize the disk image
        self.disk.normalize()
        self.disk.unit = None

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

    def make_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the map of young non-ionizing stars ...")

        # Loop over the different colour options
        #for factor in (self.config.factor_range.linear(self.config.factor_nvalues, as_list=True) + [self.config.best_factor]):
        for factor in self.config.factor_range.linear(self.config.factor_nvalues, as_list=True):

            # Calculate the non ionizing young stars map from the FUV data
            non_ionizing_stars = self.make_corrected_fuv_map(factor)

            # Add the attenuation map to the dictionary
            self.corrected_fuv_maps[factor] = non_ionizing_stars

        #best_corrected_fuv_map = self.corrected_fuv_maps[self.config.best_factor].copy()
        # Make sure all pixels of the disk-subtracted maps are larger than or equal to zero
        #best_corrected_fuv_map[best_corrected_fuv_map < 0.0] = 0.0

        # Set the best estimate of the young stars map
        #self.map = best_corrected_fuv_map

    # -----------------------------------------------------------------

    def make_corrected_fuv_map(self, factor):

        """
        This function ...
        :param factor:
        :return:
        """

        # Inform the user
        log.info("Subtracting the old stellar contribution from the map of the FUV emission with a factor of " + str(factor) + "...")

        ## Subtract old stellar contribution from FUV and MIPS 24 emission

        # From the FUV and 24 micron maps we must subtract the diffuse radiation (old stellar contribution),
        # for this we typically use an exponential disk
        # (scale length determined by GALFIT)

        flux_fuv = self.fuv.sum()

        # typisch 20% en 35% respectievelijk

        total_contribution = factor * flux_fuv

        # Subtract the disk contribution to the FUV image
        new_fuv = self.fuv - total_contribution * self.masked_disk_frame

        # Make sure all pixels of the disk-subtracted maps are larger than or equal to zero
        #new_fuv[new_fuv < 0.0] = 0.0

        # Set zero where low signal-to-noise ratio
        # new_fuv[self.fuv < self.config.non_ionizing_stars.fuv_snr_level*self.fuv_errors] = 0.0

        # Return the new FUV frame
        return new_fuv

    # -----------------------------------------------------------------

    def create_distribution_region(self):

        """
        This function ...
        :return:
        """

        disk_ellipse = self.disk_ellipse.to_pixel(self.fuv.wcs)
        inner_ellipse = disk_ellipse * self.config.histograms_annulus_range.min
        outer_ellipse = disk_ellipse * self.config.histograms_annulus_range.max
        composite = PixelCompositeRegion(outer_ellipse, inner_ellipse)
        region = PixelRegionList()
        region.append(composite)
        self.distribution_region = region

    # -----------------------------------------------------------------

    def make_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making distributions of the pixel values of the corrected FUV maps ...")

        # Create mask
        mask = self.distribution_region.to_mask(self.map.xsize, self.map.ysize)

        # Loop over the different maps
        for factor in self.corrected_fuv_maps:

            # Get the values
            values = self.corrected_fuv_maps[factor][mask]

            # Make a distribution of the pixel values indicated by the mask
            distribution = Distribution.from_values(values, bins=self.config.histograms_nbins)

            # Add the distribution to the dictionary
            self.corrected_fuv_distributions[factor] = distribution

    # -----------------------------------------------------------------

    def normalize_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Normalizing the young stellar map ...")

        # Normalize the dust map
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
