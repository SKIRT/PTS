#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.oldstars.disk Contains the DiskOldStellarMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....core.tools.logging import log
from ....core.tools import filesystem as fs
from ....magic.core.image import Image
from ....core.basics.configurable import Configurable

# -----------------------------------------------------------------

def make_map(frame, bulge):

    """
    This function ...
    :param frame:
    :param bulge:
    :return: 
    """

    # Create the maker
    maker = DiskOldStellarMapMaker()

    frames = {"standard": frame}
    bulges = {"standard": bulge}

    # Run the maker
    maker.run(frames=frames, bulges=bulges)

    # Return the map
    return maker.single_map

# -----------------------------------------------------------------

class DiskOldStellarMapMaker(Configurable):

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
        super(DiskOldStellarMapMaker, self).__init__(config, interactive)

        # -- Attributes --

        # The IRAC I1 frame in Jy
        #self.i1_jy = None

        # The IRAC I1 frame with the bulge subtracted
        #self.i1_jy_minus_bulge = None

        # The image of significance masks
        #self.significance = Image()

        # The cutoff mask
        #self.cutoff_mask = None

        self.frames = dict()
        self.bulges = dict()

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

        # 2. Load the necessary frames
        #self.load_frames()

        # 3. Calculate the significance masks
        #self.calculate_significance()

        # 4. Make the map of old stars
        self.make_maps()

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
        :return:
        """

        # Call the setup function of the base class
        super(DiskOldStellarMapMaker, self).setup(**kwargs)

        self.frames = kwargs.pop("frames")
        self.bulges = kwargs.pop("bulges")

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the maps of old stars ...")

        # Loop over the frames
        for name in self.frames:

            # Old stars = IRAC3.6 - bulge
            # From the IRAC 3.6 micron map, we must subtract the bulge component to only retain the disk emission

            # The relative contribution of the bulge to the 3.6mu emission
            #bulge_rel_contribution = self.parameters.bulge.f

            # Total flux of the IRAC 3.6mu image
            #total_flux = np.sum(self.images["3.6mu"].frames.primary)

            # Calculate factor
            #factor = bulge_rel_contribution * total_flux / np.sum(self.bulge)

            # Create the old stars map
            #old_stars = self.images["3.6mu"].frames.primary - factor * self.bulge

            #assert str(self.masked_bulge_frame.unit) == "Jy"

            frame = self.frames[name]
            bulge = self.bulges[name]

            # Subtract bulge from the IRAC I1 image
            minus_bulge = frame - bulge

            #bulge_residual = self.images["3.6mu"].frames.primary - self.disk
            #bulge_residual_path = fs.join(self.maps_intermediate_path, "bulge_residual.fits")
            #bulge_residual.save(bulge_residual_path)

            # Set the old stars map zero for pixels with low signal-to-noise in the 3.6 micron image
            #old_stars[self.irac < self.config.old_stars.irac_snr_level*self.irac_errors] = 0.0

            # Create copy
            #map = self.i1_jy_minus_bulge.copy()

            # Make sure all pixel values are larger than or equal to zero
            minus_bulge[minus_bulge < 0.0] = 0.0

            # Add
            self.maps[name] = minus_bulge

            # Mask pixels outside of the low signal-to-noise contour
            #old_stars[self.mask] = 0.0

    # -----------------------------------------------------------------

    @property
    def single_map(self):

        """
        This function ...
        :return: 
        """

        if len(self.maps) != 1: raise ValueError("Not a single map")
        return self.maps[self.maps.keys()[0]]

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
