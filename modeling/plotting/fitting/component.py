#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.plotting.fitting.chisquared Contains the ChiSquaredPlotter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..component import PlottingComponent
from ...fitting.component import FittingComponent

# -----------------------------------------------------------------

class FittingPlottingComponent(PlottingComponent, FittingComponent):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        This function ...
        """

        # Call the constructor of the base classes
        PlottingComponent.__init__(self, config)
        FittingComponent.__init__(self, config)

        self.plot_fitting_chisquared_path = None
        self.plot_fitting_distributions_path = None
        self.plot_fitting_wavelength_grids_path = None
        self.plot_fitting_dust_grids_path = None
        self.plot_fitting_celldistributions_path = None
        self.plot_fitting_runtimes_path = None
        self.plot_fitting_seds_path = None
        self.plot_fitting_images_path = None
        self.plot_fitting_geometries_path = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base classes
        PlottingComponent.setup(**kwargs)
        FittingComponent.setup(**kwargs)

        # Directory for plotting the wavelength grids
        self.plot_fitting_wavelength_grids_path = fs.create_directory_in(self.plot_fitting_path, "wavelength grids")

        # Directory for plotting the dust grids
        self.plot_fitting_dust_grids_path = fs.create_directory_in(self.plot_fitting_path, "dust grids")

        #
        self.plot_fitting_generations_path = fs.create_directory_in(self.plot_fitting_path, "generations")

        # Directory for plotting probability distributions
        self.plot_fitting_distributions_path = fs.create_directory_in(self.plot_fitting_path, "distributions")



# -----------------------------------------------------------------
