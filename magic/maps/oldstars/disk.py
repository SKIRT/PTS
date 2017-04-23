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

def make_map():

    """
    This function ...
    :return: 
    """

    maker = DiskOldStellarMapMaker()

    maker.run()

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
        self.i1_jy = None

        # The IRAC I1 frame with the bulge subtracted
        self.i1_jy_minus_bulge = None

        # The map of the old stars
        self.map = None

        # The image of significance masks
        self.significance = Image()

        # The cutoff mask
        self.cutoff_mask = None

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
        if self.config.write: self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(DiskOldStellarMapMaker, self).setup(**kwargs)

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

        assert str(self.masked_bulge_frame.unit) == "Jy"

        # Subtract bulge from the IRAC I1 image
        self.i1_jy_minus_bulge = self.i1_jy - self.masked_bulge_frame

        #bulge_residual = self.images["3.6mu"].frames.primary - self.disk
        #bulge_residual_path = fs.join(self.maps_intermediate_path, "bulge_residual.fits")
        #bulge_residual.save(bulge_residual_path)

        # Set the old stars map zero for pixels with low signal-to-noise in the 3.6 micron image
        #old_stars[self.irac < self.config.old_stars.irac_snr_level*self.irac_errors] = 0.0

        # Create copy
        self.map = self.i1_jy_minus_bulge.copy()

        # Make sure all pixel values are larger than or equal to zero
        self.map[self.map < 0.0] = 0.0

        # Mask pixels outside of the low signal-to-noise contour
        #old_stars[self.mask] = 0.0

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

        # Write the IRAC I1 image with the bulge subtracted
        self.write_i1_minus_bulge()

        # Write the map of old stars
        self.write_map()

        # Write the significance mask
        self.write_significance_masks()

    # -----------------------------------------------------------------

    def write_i1_minus_bulge(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the IRAC I1 image with the bulge subtracted ...")

        # Determine path
        path = fs.join(self.maps_old_path, "IRAC I1 min bulge.fits")

        # Write
        self.i1_jy_minus_bulge.saveto(path)

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
