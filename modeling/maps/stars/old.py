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
from ..component import MapsComponent
from ....core.tools import filesystem as fs

# -----------------------------------------------------------------

class OldStellarMapMaker(MapsComponent):

    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(OldStellarMapMaker, self).__init__(config)

        # -- Attributes --

        # The IRAC I1 frame in Jy
        self.i1_jy = None

        # The IRAC I1 frame with the bulge subtracted
        self.i1_jy_minus_bulge = None

        # The map of the old stars
        self.map = None

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

        # 3. Make the map of old stars
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
        super(OldStellarMapMaker, self).setup()

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

        assert str(self.bulge_frame.unit) == "Jy"

        # Subtract bulge from the IRAC I1 image
        self.i1_jy_minus_bulge = self.i1_jy - self.bulge_frame

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
        self.i1_jy_minus_bulge.save(path)

    # -----------------------------------------------------------------

    def write_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the map of old stars ...")

        # Write
        self.map.save(self.old_stellar_map_path)

# -----------------------------------------------------------------
