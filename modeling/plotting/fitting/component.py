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
from ....core.tools import filesystem as fs

# -----------------------------------------------------------------

class FittingPlottingComponent(PlottingComponent, FittingComponent):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param interactive:
        """

        # Call the constructor of the base classes
        PlottingComponent.__init__(self, no_config=True)
        FittingComponent.__init__(self, *args, **kwargs)

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
        PlottingComponent.setup(self, **kwargs)
        FittingComponent.setup(self, **kwargs)

        # Make paths
        self.plot_fitting_chisquared_path = fs.create_directory_in(self.plot_fitting_path, "chi squared")
        self.plot_fitting_distributions_path = fs.create_directory_in(self.plot_fitting_path, "distributions")
        self.plot_fitting_wavelength_grids_path = fs.create_directory_in(self.plot_fitting_path, "wavelength grids")
        self.plot_fitting_dust_grids_path = fs.create_directory_in(self.plot_fitting_path, "dust grids")
        self.plot_fitting_celldistributions_path = fs.create_directory_in(self.plot_fitting_path, "dust cell distributions")
        self.plot_fitting_runtimes_path = fs.create_directory_in(self.plot_fitting_path, "runtimes")
        self.plot_fitting_seds_path = fs.create_directory_in(self.plot_fitting_path, "seds")
        self.plot_fitting_images_path = fs.create_directory_in(self.plot_fitting_path, "images")
        self.plot_fitting_geometries_path = fs.create_directory_in(self.plot_fitting_path, "geometries")

# -----------------------------------------------------------------
