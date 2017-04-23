#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.maps.youngstars.young Contains the YoungStellarMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....core.tools.logging import log
from ....core.tools import filesystem as fs
from ....core.basics.distribution import Distribution
from ....magic.region.composite import PixelCompositeRegion
from ....magic.region.list import PixelRegionList
from ....core.basics.configurable import Configurable

# -----------------------------------------------------------------

def make_map(fuv, fuv_atttenuation, old, factor):

    """
    This function ...
    :return: 
    """

    # Create the map maker
    maker = YoungStellarMapsMaker()

    # Set input
    factors = [factor]
    fuv_attenuations = {"standard": fuv_atttenuation}

    # Run the map maker
    maker.run(fuv=fuv, fuv_attenuations=fuv_attenuations, old=old, factors=factors)

    # Return the map
    return maker.single_map

# -----------------------------------------------------------------

def make_maps(fuv, fuv_attenuation, old, factors):

    """
    THis function ...
    :param fuv: 
    :param fuv_attenuation: 
    :param old: 
    :param factors: 
    :return: 
    """

    # Create the map maker
    maker = YoungStellarMapsMaker()

    # Set input
    fuv_attenuations = {"standard": fuv_attenuation}

    # Run the map maker
    maker.run(fuv=fuv, fuv_attenuations=fuv_attenuations, old=old, factors=factors)

    # Return the maps
    return maker.maps

# -----------------------------------------------------------------

class YoungStellarMapsMaker(Configurable):

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
        super(YoungStellarMapsMaker, self).__init__(config, interactive)

        # -- Attributes --

        # The input FUV and FUV error maps
        self.fuv = None
        self.fuv_errors = None

        self.old = None

        self.fuv_attenuations = None

        self.factors = None

        # The maps of the corrected FUV emission
        #self.corrected_fuv_maps = dict()

        # The distributions of corrected FUV pixel values
        #self.corrected_fuv_distributions = dict()

        # The map of young stars
        #self.map = None

        # The path to the maps/young/fuv directory
        #self.maps_young_fuv_path = None

        # Region of area taken for calculating distribution of pixel values
        #self.distribution_region = None

        # The image of significance masks
        #self.significance = Image()

        # The cutoff mask
        #self.cutoff_mask = None

        # The transparent FUV maps
        self.transparent = dict()

        # The maps
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

        # 3. Calculate the significance masks
        #self.calculate_significance()

        # 3. Make the map of young stars
        self.make_maps()

        # ...
        #self.create_distribution_region()
        #self.make_distributions()

        # 4. Normalize the map
        #self.normalize_map()

        # Create the cutoff mask
        #self.make_cutoff_mask()

        # 5. Cut-off map
        #self.cutoff_map()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(YoungStellarMapsMaker, self).setup(**kwargs)

        # Get input
        self.fuv = kwargs.pop("fuv")
        self.fuv_errors = kwargs.pop("fuv_errors", None)
        self.old = kwargs.pop("old")
        self.fuv_attenuations = kwargs.pop("fuv_attenuations")

        # Set factors
        self.factors = kwargs.pop("factors")

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the map of young non-ionizing stars ...")

        # Correct for internal attenuation
        self.correct_for_attenuation()

        # Subtract the contribution of old stars
        self.subtract_old_contribution()

    # -----------------------------------------------------------------

    def correct_for_attenuation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Correcting for attenuation ...")

        # attenuation = -2.5 log (total / transparent)

        # Loop over the different attenuation maps
        for name in self.fuv_attenuations:

            # Get the attenuation map
            attenuation = self.fuv_attenuations[name]

            # Calculate the transparent FUV flux
            exponent = attenuation / 2.5
            transparent = self.fuv * 10**exponent

            # Set the corrected map
            self.transparent[name] = transparent

    # -----------------------------------------------------------------

    def subtract_old_contribution(self):

        """
        Thisnfunction ...
        :return:
        """

        # Inform the user
        log.info("Subtracting the old contribution ...")

        # Loop over the transparent maps
        for name in self.transparent:

            # Get the transparent map
            transparent = self.transparent[name]

            # Loop over the different factors
            for factor in self.factors:

                # Calculate the non ionizing young stars map from the FUV data
                young_stars = make_corrected_fuv_map(transparent, self.old, factor)

                # Determine name
                key = name + "_" + repr(factor)

                # Add the attenuation map to the dictionary
                self.maps[key] = young_stars

    # -----------------------------------------------------------------

    @property
    def single_map(self):

        """
        This fucntion ...
        :return: 
        """

        if len(self.maps) != 1: raise ValueError("Not a single map")
        return self.maps[self.maps.keys()[0]]

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

def make_corrected_fuv_map(fuv, old, factor):

    """
    This function ...
    :param fuv:
    :param old:
    :param factor:
    :return:
    """

    # Inform the user
    log.info("Subtracting the old stellar contribution from the map of the FUV emission with a factor of " + str(factor) + "...")

    ## Subtract old stellar contribution from FUV and MIPS 24 emission

    # From the FUV and 24 micron maps we must subtract the diffuse radiation (old stellar contribution),
    # for this we typically use an exponential disk
    # (scale length determined by GALFIT)

    flux_fuv = fuv.sum()

    # typisch 20% en 35% respectievelijk

    total_contribution = factor * flux_fuv

    # Subtract the disk contribution to the FUV image
    new_fuv = fuv - total_contribution * old

    # Make sure all pixels of the disk-subtracted maps are larger than or equal to zero
    #new_fuv[new_fuv < 0.0] = 0.0

    # Set zero where low signal-to-noise ratio
    # new_fuv[self.fuv < self.config.non_ionizing_stars.fuv_snr_level*self.fuv_errors] = 0.0

    # Return the new FUV frame
    return new_fuv

# -----------------------------------------------------------------
