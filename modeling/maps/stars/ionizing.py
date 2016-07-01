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

# Import the relevant PTS classes and modules
from ....core.tools.logging import log

# -----------------------------------------------------------------

class IonizingStellarMapMaker(object):

    """
    This class...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(IonizingStellarMapMaker, self).__init__()

        # -- Attributes --

        self.map = None

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # ...

        #
        self.make_map()

        # Return the map
        return self.map

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def make_map(self):

        """
        This function ...
        :return:
        """

        # H-ALPHA HAS BEEN CONVERTED TO LSUN (ABOVE)

        # Young ionizing stars = Ha + 0.031 x MIPS24

        # CALCULATE THE IONIZING STARS MAP BASED ON THE CONVERTED H ALPHA AND THE DISK-SUBTRACTED 24 MICRON IMAGE

        # Calculate the young stellar contribution to the 24 micron image
        mips_young_stars = self.get_mips_young_stars_map()

        # Save the mips_young_stars map
        mips_young_path = fs.join(self.maps_intermediate_path, "24mu_young.fits")
        mips_young_stars.save(mips_young_path)

        # Calculate ionizing stars map and ratio
        ionizing = self.images["Halpha"].frames.primary + 0.031 * mips_young_stars

        # ionizing_ratio = self.ha / (0.031*mips_young_stars)

        # MASK NEGATIVE AND LOW SIGNAL-TO-NOISE PIXELS

        # Set pixels to zero with low signal-to-noise in the H Alpha image
        # ionizing[self.ha < self.config.ionizing_stars.ha_snr_level*self.ha_errors] = 0.0
        # ionizing_ratio[self.ha < self.config.ionizing_stars.ha_snr_level*self.ha_errors] = 0.0

        # Set pixels to zero with low signal-to-noise in the 24 micron image
        # ionizing[self.mips < self.config.ionizing_stars.mips_snr_level*self.mips_errors] = 0.0
        # ionizing_ratio[self.mips < self.config.ionizing_stars.mips_snr_level*self.mips_errors] = 0.0

        # Make sure all pixel values are larger than or equal to zero
        # ionizing[ionizing < 0.0] = 0.0
        # ionizing_ratio[ionizing < 0.0] = 0.0

        # New
        ionizing[self.cutoff_masks["Halpha"]] = 0.0

        # Set the ionizing stars map
        self.map = ionizing

    # -----------------------------------------------------------------

    def get_mips_young_stars_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Subtracting the old stellar contribution from the 24 micron emission ...")

        ## Subtract old stellar contribution from FUV and MIPS 24 emission

        #     From the FUV and 24 micron maps we must subtract the diffuse radiation (old stellar contribution),
        #     for this we typically use an exponential disk
        #     (scale length detemermined by GALFIT)

        ## MIPS HAS BEEN CONVERTED TO LSUN (ABOVE)

        flux_mips = np.sum(self.images["24mu"].frames.primary)

        # typisch 20% en 35% respectievelijk
        # 48% voor MIPS 24 komt van Lu et al. 2014

        factor = 0.48 * flux_mips / np.sum(self.disk)

        # Subtract the disk contribution to the 24 micron image
        new_mips = self.images["24mu"].frames.primary - factor * self.disk

        # Make sure all pixels of the disk-subtracted maps are larger than or equal to zero
        new_mips[new_mips < 0.0] = 0.0

        # Set zero where low signal-to-noise ratio
        # new_mips[self.mips < self.config.ionizing_stars.mips_young_stars.mips_snr_level*self.mips_errors] = 0.0

        # Return the new 24 micron frame
        return new_mips

# -----------------------------------------------------------------
