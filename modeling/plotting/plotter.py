#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.plotting.plotter Contains the Plotter class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import PlottingComponent
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
from .data import DataPlotter
from .preparation import PreparationPlotter
from .decomposition import DecompositionPlotter
from .truncation import TruncationPlotter
from .photometry import PhotometryPlotter
from .fitting import FittingPlotter
from .maps import MapsPlotter
from ...magic.plot.imagegrid import StandardImageGridPlotter, ResidualImageGridPlotter

# -----------------------------------------------------------------

class Plotter(PlottingComponent):

    """
    This class ...
    """

    # -----------------------------------------------------------------

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(Plotter, self).__init__(config)

        # -- Attributes --

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new Plotter instance
        plotter = cls(arguments.config)

        # Set the modeling path
        plotter.config.path = arguments.path

        # Set the modeling step for which to make the plots
        plotter.config.step = arguments.step

        # Return the new instance
        return plotter

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Make a plot of the data
        if self.config.step == "data": self.plot_data()

        # 3. Make a plot of the preparation step
        elif self.config.step == "preparation": self.plot_preparation()

        # 4. Make a plot of the decomposition step
        elif self.config.step == "decomposition": self.plot_decomposition()

        # 5. Make a plot of the truncation step
        elif self.config.step == "truncation": self.plot_truncation()

        # 6. Make a plot of the photometry step
        elif self.config.step == "photometry": self.plot_photometry()

        # 7. Make a plot of the map making step
        elif self.config.step == "maps": self.plot_maps()

        # 8. Make a plot of the fitting step
        elif self.config.step == "fit": self.plot_fit()

        # Invalid modelling step
        else: raise ValueError("Invalid modelling step")

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(Plotter, self).setup()

    # -----------------------------------------------------------------

    def plot_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        plotter = DataPlotter()
        plotter.run()

    # -----------------------------------------------------------------

    def plot_preparation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        plotter = PreparationPlotter()
        plotter.run()

    # -----------------------------------------------------------------

    def plot_decomposition(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        plotter = DecompositionPlotter()
        plotter.run()

    # -----------------------------------------------------------------

    def plot_truncation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        plotter = TruncationPlotter()
        plotter.run()

    # -----------------------------------------------------------------

    def plot_photometry(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        plotter = PhotometryPlotter()
        plotter.run()

    # -----------------------------------------------------------------

    def plot_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        plotter = MapsPlotter()
        plotter.run()

    # -----------------------------------------------------------------

    def plot_fit(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        plotter = FittingPlotter()
        plotter.run()

# -----------------------------------------------------------------
