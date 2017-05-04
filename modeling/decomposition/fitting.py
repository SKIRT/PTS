#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.decomposition.fitting Contains the FittingDecomposer class.

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

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(FittingDecomposer, self).__init__(*args, **kwargs)

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

# From https://groups.google.com/forum/#!topic/astropy-dev/W0GzYoSvjF4
def moffat_fitting():

    from astropy.modeling.models import custom_model
    # Define model
    @custom_model
    def Elliptical_Moffat2D(x, y, N_sky = 1., amplitude = 1., phi=0., power = 1., x_0 = 0., y_0 = 0., width_x = 1., width_y = 1.):

        c = np.cos(phi)
        s = np.sin(phi)
        A = (c / width_x) ** 2 + (s / width_y)**2
        B = (s / width_x) ** 2 + (c/ width_y)**2
        C = 2 * s * c * (1/width_x**2 - 1/width_y**2)
        denom = (1 + A * (x-x_0)**2 + B * (y-y_0)**2 + C*(x-x_0)*(y-y_0))**power
        return N_sky + amplitude / denom

    moffat_init  = Elliptical_Moffat2D(amplitude=peak_targ, x_0=cbox_size, y_0=cbox_size)
    fit_y, fit_x = np.mgrid[:2*cbox_size, :2*cbox_size]
    fit_moffat   = fitting.LevMarLSQFitter()
    moffat       = fit_moffat(moffat_init, fit_x, fit_y, image_targ)
    #moffat

# -----------------------------------------------------------------
