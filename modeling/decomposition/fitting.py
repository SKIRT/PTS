#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.decomposition TO DO

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.modeling.models import Sersic2D

# Import the relevant PTS classes and modules
from .component import DecompositionComponent

# -----------------------------------------------------------------

class FittingDecompositionParameters(DecompositionComponent):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(FittingDecompositionParameters, self).__init__(config)

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        sersic_model = Sersic2D(amplitude=sersic_amplitide, r_eff=sersic_r_eff, n=sersic_n, x_0=sersic_x_0, y_0=sersic_y_0, ellip=sersic_ellip, theta=sersic_theta)
        sersic_map = sersic_model(sersic_x, sersic_y)

# -----------------------------------------------------------------
