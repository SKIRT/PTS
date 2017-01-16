#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.basics.pixelscale Contains the Pixelscale class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
from astropy.units import Unit

# Import standard modules
from .vector import Extent

# -----------------------------------------------------------------

class Pixelscale(Extent):

    """
    This class ...
    """

    def __init__(self, x, y=None):

        """
        The constructor ...
        """

        # Remove '/pix' if present
        if y is None: y = x

        # Convert to just arcsec or just
        x = only_angle(x)
        y = only_angle(y)

        # Call the constructor of the base class
        super(Pixelscale, self).__init__(x, y)

    # -----------------------------------------------------------------

    @property
    def abs_x(self):

        """
        This function ...
        :return:
        """

        return abs(self.x)

    # -----------------------------------------------------------------

    @property
    def abs_y(self):

        """
        This function ...
        :return:
        """

        return abs(self.y)

    # -----------------------------------------------------------------

    @property
    def average(self):

        """
        This function ...
        :return:
        """

        x_pixelscale = abs(self.x.to("arcsec"))
        y_pixelscale = abs(self.y.to("arcsec"))

        # Check if the x pixelscale and y pixelscale are close
        if not np.isclose(x_pixelscale.value, y_pixelscale.value, rtol=0.001):

            # Warn about the difference in x and y pixelscale
            from ...core.tools.logging import log
            log.warning("Averaging the pixelscale over the x and y direction may not be a good approximation:")
            log.warning("  * x pixelscale (absolute value) = " + str(x_pixelscale))
            log.warning("  * y pixelscale (absolute value) = " + str(y_pixelscale))

        # Return a single value for the pixelscale in arcseconds
        return 0.5 * (x_pixelscale + y_pixelscale)

        # Other way from wcs:
        #np.mean(np.abs(np.diagonal(img_wcs.pixel_scale_matrix)))

    # -----------------------------------------------------------------

    @property
    def solid_angle(self):

        """
        This function ...
        :return:
        """

        solid_angle = (self.abs_x * self.abs_y).to("sr")
        return solid_angle

# -----------------------------------------------------------------

def only_angle(quantity):

    """
    This function ...
    :return:
    """

    if "pix" in quantity.unit.bases: return quantity * Unit("pix")
    elif len(quantity.unit.bases) == 1: return quantity
    else: raise ValueError("Don't know what to do with " + str(quantity) + " to convert to angular pixelscale")

# -----------------------------------------------------------------
