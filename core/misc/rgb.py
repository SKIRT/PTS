#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.misc.rgb Contains the RGBImageMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
from astropy.visualization import make_lupton_rgb
from reproject import reproject_interp
from astropy.io import fits

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable

# -----------------------------------------------------------------

def make_example_rgbs(self):

    """
    This function ...
    :return:
    """

    # Read in the three images downloaded from here:
    # g: http://dr13.sdss.org/sas/dr13/eboss/photoObj/frames/301/1737/5/frame-g-001737-5-0039.fits.bz2
    # r: http://dr13.sdss.org/sas/dr13/eboss/photoObj/frames/301/1737/5/frame-r-001737-5-0039.fits.bz2
    # i: http://dr13.sdss.org/sas/dr13/eboss/photoObj/frames/301/1737/5/frame-i-001737-5-0039.fits.bz2
    g = fits.open('frame-g-001737-5-0039.fits.bz2')[0]
    r = fits.open('frame-r-001737-5-0039.fits.bz2')[0]
    i = fits.open('frame-i-001737-5-0039.fits.bz2')[0]

    # remap r and i onto g
    r_new, r_mask = reproject_interp(r, g.header)
    i_new, i_mask = reproject_interp(i, g.header)

    # zero out the unmapped values
    i_new[np.logical_not(i_mask)] = 0
    r_new[np.logical_not(r_mask)] = 0

    # red=i, green=r, blue=g
    # make a file with the default scaling
    rgb_default = make_lupton_rgb(i_new, r_new, g.data, filename="ngc6976-default.jpeg")
    # this scaling is very similar to the one used in Lupton et al. (2004)
    rgb = make_lupton_rgb(i_new, r_new, g.data, Q=10, stretch=0.5, filename="ngc6976.jpeg")

# -----------------------------------------------------------------

class RGBImageMaker(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(RGBImageMaker, self).__init__(*args, **kwargs)

        # -- Attributes --



    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(RGBImageMaker, self).setup(**kwargs)

# -----------------------------------------------------------------
