#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.plotting.fitting Contains the FittingPlotter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import PlottingComponent
from ..fitting.component import FittingComponent

# Import the relevant
from ...component.component import load_modeling_configuration
from .chisquared import ChiSquardedPLotter
from .distributions import DistributionsPLotter
from .wavelengthgrids import WavelengthGridsPlotter
from .dustgrids import DustGridsPlotter
from .celldistributions import CellDistributionsPlotter
from .runtimes import RuntimesPlotter
from .seds import SEDsPLotter
from .images import ImagesPlotter
from .geometries import GeometriesPlotter
from ....core.basics.log import log
from ....core.basics.configurable import Configurable

# -----------------------------------------------------------------

def get_features(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    # Load the modeling configuration
    modeling_configuration = load_modeling_configuration(modeling_path)

    features = dict()
    features["chi-squared"] = "chi-squared"
    features["distributions"] = "distributions"
    features["wavelength grids"] = "wavelength grids"

    if modeling_configuration.type == "galaxy": features["dust grids"] = "dust grids"
    if modeling_configuration.type == "galaxy": features["cell distribution"] = "cell distribution"

    features["runtimes"] = "runtimes"
    features["seds"] = "seds"
    features["images"] = "images"

    if modeling_configuration.type == "galaxy": features["geometries"] = "geometries"

    # Return
    return features

# -----------------------------------------------------------------

class FittingPlotter(Configurable):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base classes
        #PlottingComponent.__init__(self, config)
        #FittingComponent.__init__(self, config)

        super(FittingPlotter, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Plot chi-squared
        if "chi-squared" in self.config.features: self.plot_chi_squared()

        # Distributions
        if "distribution" in self.config.features: self.plot_distributions()

        # Wavelength grids
        if "wavelength grids" in self.config.features: self.plot_wavelength_grids()

        # Dust grids
        if "dust grids" in self.config.features: self.plot_dust_grids()

        # Cell distributions
        if "cell distribution" in self.config.features: self.plot_dust_cell_distributions()

        # Runtimes
        if "runtimes" in self.config.features: self.plot_runtimes()

        # SEDs
        if "seds" in self.config.features: self.plot_seds()

        # Images
        if "images" in self.config.features: self.plot_images()

        # Geometries
        if "geometries" in self.config.features: self.plot_geometries()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(FittingPlotter, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def plot_chi_squared(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting chi squared ...")

        plotter = ChiSquardedPLotter()
        plotter.config.path = self.config.path
        plotter.run()

    # -----------------------------------------------------------------

    def plot_distributions(self):

        """
        THis function ...
        :return:
        """

        # Inform the user
        log.info("Plotting distributions ...")

        plotter = DistributionsPLotter()
        plotter.config.path = self.config.path
        plotter.run()

    # -----------------------------------------------------------------

    def plot_wavelength_grids(self):

        """
        This fucntion ...
        :return:
        """

        # Inform the user
        log.info("Plotting wavelength grids ...")

        plotter = WavelengthGridsPlotter()
        plotter.config.path = self.config.path
        plotter.run()

    # -----------------------------------------------------------------

    def plot_dust_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting dust grids ...")

        # Create plotter
        plotter = DustGridsPlotter()
        plotter.config.path = self.config.path
        plotter.run()

    # -----------------------------------------------------------------

    def plot_dust_cell_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting dust cell distributions ...")

        # Create plotter
        plotter = CellDistributionsPlotter()
        plotter.config.path = self.config.path
        plotter.run()

    # -----------------------------------------------------------------

    def plot_runtimes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the runtimes ...")

        # Create plotter
        plotter = RuntimesPlotter()
        plotter.config.path = self.config.path
        plotter.run()

    # -----------------------------------------------------------------

    def plot_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting SEDs ...")

        plotter = SEDsPLotter()
        plotter.config.path = self.config.path
        plotter.run()

    # -----------------------------------------------------------------

    def plot_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting images ...")

        plotter = ImagesPlotter()
        plotter.config.path = self.config.path
        plotter.run()

    # -----------------------------------------------------------------

    def plot_geometries(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting geometries ...")

        plotter = GeometriesPlotter()
        plotter.config.path = self.config.path
        plotter.run()

# -----------------------------------------------------------------
