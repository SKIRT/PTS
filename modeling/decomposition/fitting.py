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

class FittingDecomposer(DecompositionComponent):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(FittingDecomposer, self).__init__(config)

        # The dictionary of components
        self.components = dict()

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the input images
        self.load_images()

        # 3. Do the fitting
        self.fit()

        # 4. Create the component models
        self.create_components()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(FittingDecomposer, self).setup()

    # -----------------------------------------------------------------

    def load_images(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def fit(self):

        """
        This function ...
        :return:
        """

        sersic_model = Sersic2D(amplitude=sersic_amplitide, r_eff=sersic_r_eff, n=sersic_n, x_0=sersic_x_0, y_0=sersic_y_0, ellip=sersic_ellip, theta=sersic_theta)
        sersic_map = sersic_model(sersic_x, sersic_y)

    # -----------------------------------------------------------------

    def create_components(self):

        """
        This function ...
        :return:
        """

# -----------------------------------------------------------------
