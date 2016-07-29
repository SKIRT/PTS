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

# -----------------------------------------------------------------

class YoungStellarMapMaker(MapsComponent):

    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(YoungStellarMapMaker, self).__init__(config)

        # -- Attributes --

        # The input FUV and FUV error maps
        self.fuv = None
        self.fuv_errors = None

        # The NORMALIZED (to unity) DISK IMAGE
        self.disk = None

        # The maps of the corrected FUV emission
        self.corrected_fuv_maps = dict()

        # The map of young stars
        self.map = None

        # The path to the maps/young/fuv directory
        self.maps_young_fuv_path = None

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

        # 3. Make the map of young stars
        self.make_map()

        # 4. Normalize the map
        self.normalize_map()

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
        self.fuv_errors = self.dataset.get_errors("GALEX FUV") # in original MJy/sr units

    # -----------------------------------------------------------------

    def load_disk(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the disk image ...")

        # Get disk frame
        self.disk = self.disk_frame

        # Normalize the disk image
        self.disk.normalize()
        self.disk.unit = None

    # -----------------------------------------------------------------

    def make_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the map of young non-ionizing stars ...")

        # Loop over the different colour options
        for factor in (self.config.factor_range.linear(self.config.factor_nvalues, as_list=True) + [self.config.best_factor]):

            # Calculate the non ionizing young stars map from the FUV data
            non_ionizing_stars = self.make_corrected_fuv_map(factor)

            # Add the attenuation map to the dictionary
            self.corrected_fuv_maps[factor] = non_ionizing_stars

        # Set the best estimate of the young stars map
        self.map = self.corrected_fuv_maps[self.config.best_factor]

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
        new_fuv = self.fuv - total_contribution * self.disk_frame

        # Make sure all pixels of the disk-subtracted maps are larger than or equal to zero
        new_fuv[new_fuv < 0.0] = 0.0

        # Set zero where low signal-to-noise ratio
        # new_fuv[self.fuv < self.config.non_ionizing_stars.fuv_snr_level*self.fuv_errors] = 0.0

        # Return the new FUV frame
        return new_fuv

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

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the corrected FUV maps
        self.write_fuv_maps()

        # Write the final young stellar map
        self.write_map()

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
            self.corrected_fuv_maps[factor].save(path)

    # -----------------------------------------------------------------

    def write_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the map of young stars ...")

        # Write
        self.map.save(self.young_stellar_map_path)

# -----------------------------------------------------------------
