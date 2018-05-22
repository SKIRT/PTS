#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.plotting.component Contains the PlottingComponent class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta, abstractmethod

# Import the relevant PTS classes and modules
from ..component.component import ModelingComponent
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

class PlottingComponent(ModelingComponent):
    
    """
    This class...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(PlottingComponent, self).__init__(*args, **kwargs)

        # -- Attributes --

        # Plotting subdirectories
        self.plot_data_path = None
        self.plot_preparation_path = None
        self.plot_decomposition_path = None
        self.plot_truncation_path = None
        self.plot_photometry_path = None
        self.plot_maps_path = None
        self.plot_fitting_path = None
        self.plot_analysis_path = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(PlottingComponent, self).setup(**kwargs)

        # Set the path to the plot/data directory
        self.plot_data_path = fs.create_directory_in(self.plot_path, "data")

        # Set the path to the plot/prep directory
        self.plot_preparation_path = fs.create_directory_in(self.plot_path, "prep")

        # Set the path to the plot/components directory
        self.plot_decomposition_path = fs.create_directory_in(self.plot_path, "components")

        # Set the path to the plot/trunc directory
        self.plot_truncation_path = fs.create_directory_in(self.plot_path, "trunc")

        # Set the path to the plot/phot directory
        self.plot_photometry_path = fs.create_directory_in(self.plot_path, "phot")

        # Set the path to the plot/maps directory
        self.plot_maps_path = fs.create_directory_in(self.plot_path, "maps")

        # Set the path to the plot/fit directory
        self.plot_fitting_path = fs.create_directory_in(self.plot_path, "fit")

        # Set the path to the plot/analysis directory
        self.plot_analysis_path = fs.create_directory_in(self.plot_path, "analysis")

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Load the data to be plotted
        self.load()

        # 3. Plot
        self.plot()

    # -----------------------------------------------------------------

    @abstractmethod
    def load(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractmethod
    def plot(self):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------
