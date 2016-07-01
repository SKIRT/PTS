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

# -----------------------------------------------------------------

class YoungStellarMapMaker(object):

    """
    This class...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(YoungStellarMapMaker, self).__init__()

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

        # Calculate the non ionizing young stars map from the FUV data
        non_ionizing_stars = self.get_fuv_young_stars_map()

        # Set the young stars map
        self.map = non_ionizing_stars

    # -----------------------------------------------------------------

    def get_fuv_young_stars_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Subtracting the old stellar contribution from the FUV emission ...")

        ## Subtract old stellar contribution from FUV and MIPS 24 emission

        #     From the FUV and 24 micron maps we must subtract the diffuse radiation (old stellar contribution),
        #     for this we typically use an exponential disk
        #     (scale length detemermined by GALFIT)

        flux_fuv = np.sum(self.images["FUV"].frames.primary)

        # typisch 20% en 35% respectievelijk
        # 48% voor MIPS 24 komt van Lu et al. 2014

        factor = 0.2 * flux_fuv / np.sum(self.disk)

        # Subtract the disk contribution to the FUV image
        new_fuv = self.images["FUV"].frames.primary - factor * self.disk

        # Make sure all pixels of the disk-subtracted maps are larger than or equal to zero
        new_fuv[new_fuv < 0.0] = 0.0

        # Set zero where low signal-to-noise ratio
        # new_fuv[self.fuv < self.config.non_ionizing_stars.fuv_snr_level*self.fuv_errors] = 0.0

        # Return the new FUV frame
        return new_fuv

# -----------------------------------------------------------------
