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

        # Call the setup function of the base class
        super(OldStellarMapMaker, self).setup()

    # -----------------------------------------------------------------

    def make_map(self):

        """
        This function ...
        :return:
        """

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

        # Convert the 3.6 micron image from MJy/sr to Jy/sr
        conversion_factor = 1.0
        conversion_factor *= 1e6

        # Convert the 3.6 micron image from Jy / sr to Jy / pixel
        pixelscale = self.images["3.6mu"].average_pixelscale
        pixel_factor = (1.0/pixelscale**2).to("pix2/sr").value
        conversion_factor /= pixel_factor
        self.images["3.6mu"] *= conversion_factor
        self.images["3.6mu"].unit = "Jy"

        i1_jy_path = fs.join(self.maps_intermediate_path, "i1_jy.fits")
        self.images["3.6mu"].save(i1_jy_path)

        # Subtract bulge
        #old_stars = self.images["3.6mu"].frames.primary - (self.bulge * 1.5)
        old_stars = self.images["3.6mu"].frames.primary - self.bulge

        bulge_residual = self.images["3.6mu"].frames.primary - self.disk
        bulge_residual_path = fs.join(self.maps_intermediate_path, "bulge_residual.fits")
        bulge_residual.save(bulge_residual_path)

        # Set the old stars map zero for pixels with low signal-to-noise in the 3.6 micron image
        #old_stars[self.irac < self.config.old_stars.irac_snr_level*self.irac_errors] = 0.0

        # Make sure all pixel values are larger than or equal to zero
        old_stars[old_stars < 0.0] = 0.0

        # Mask pixels outside of the low signal-to-noise contour
        #old_stars[self.mask] = 0.0

        # Set the old stars map
        self.map = old_stars

# -----------------------------------------------------------------
