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
from .component import FittingPlottingComponent
from ....core.basics.log import log
from ....core.tools import filesystem as fs
from ....core.basics.log import log
from ....core.simulation.logfile import LogFile

# -----------------------------------------------------------------

class CellDistributionsPlotter(FittingPlottingComponent):

    """
    This function ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(CellDistributionsPlotter, self).__init__(*args, **kwargs)

        # The distribution of dust cells per level
        self.dust_cell_trees = []

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Plot
        self.plot()

    # -----------------------------------------------------------------

    def load_dust_cell_trees(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the dust cell tree information ...")

        # Loop over the directories inside the fit/dust grids directory
        for path in fs.directories_in_path(self.fitting_run.dust_grids_path, sort=int):

            # Determine the path to the log file of the dust grid generating simulation
            log_file_path = fs.join(path, self.object_name + "_log.txt")

            # If the log file does not exist, skip
            if not fs.is_file(log_file_path): continue

            # Open the log file
            log_file = LogFile(log_file_path)

            # Get the distribution of cells per level of the tree
            self.dust_cell_trees.append(log_file.tree_leaf_distribution)

    # -----------------------------------------------------------------

    def plot_dust_cell_distribution(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the dust cell distribution ...")

        # Loop over the distributions
        index = 0
        for distribution in self.dust_cell_trees:
            title = "Dust cells in each tree level"
            path = fs.join(self.plot_fitting_dust_grids_path, "cells_tree_" + str(index) + ".pdf")

            # Plot
            distribution.plot(title=title, path=path)

            # Increment the index
            index += 1

# -----------------------------------------------------------------
