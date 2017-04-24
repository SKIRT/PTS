#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.strategy Determine the modelling strategy for a certain galaxy.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

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

    def calculate_significance(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the significance masks ...")

        # Get the significance mask
        if self.config.fuv_significance > 0: self.significance.add_mask(self.dataset.get_significance_mask("GALEX FUV", self.config.fuv_significance), "GALEX_FUV")
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
