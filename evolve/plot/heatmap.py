#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.evolve.plot.heatmap Contains the HeatMapPlotter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import matplotlib.pyplot as plt

# Import the relevant PTS classes and modules
from .plotter import Plotter
from ...core.basics.log import log
from ..analyse.database import get_generations, get_scores, get_fitnesses

# -----------------------------------------------------------------

class HeatMapPlotter(Plotter):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(HeatMapPlotter, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Get the data
        self.get_data()

        # Plot
        self.plot()

        # Write
        self.write()

    # -----------------------------------------------------------------

    def get_data(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Getting the data ...")

        self.data = dict()

        # Loop over the runs
        for run_id in self.runs:

            data = []

            # Loop over the generations
            for generation in get_generations(self.database, run_id):

                # Get the scores
                if self.config.fitness: scores = get_fitnesses(self.database, run_id, generation)
                else: scores = get_scores(self.database, run_id, generation)

                data.append(scores)

            # Set the data for this run
            self.data[run_id] = data

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Plotting ...")

        # Loop over the runs
        for run_id in self.data: self.plot_heatmap(run_id)

    # -----------------------------------------------------------------

    def plot_heatmap(self, run_id):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Plotting ...")

        plt.figure()

        # Get the data
        data = self.data[run_id]

        # Create
        plt.imshow(data, aspect="auto", interpolation="gaussian", cmap="viridis")

        self.set_title("Plot of population scores along the generations for run '" + run_id + "'")

        plt.xlabel('Population')
        plt.ylabel('Generations')

        plt.grid(True)
        plt.colorbar()

        # Show
        self.write_or_show("heatmap")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Writing ...")

        # Write the data
        self.write_data()

# -----------------------------------------------------------------
