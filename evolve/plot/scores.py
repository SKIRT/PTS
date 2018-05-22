#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.evolve.plot.scores Contains the ScoresPlotter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import matplotlib.pyplot as plt

# Import the relevant PTS classes and modules
from .plotter import Plotter
from ...core.basics.log import log
from ..analyse.database import get_generations, get_statistics
from ...core.basics.map import Map

# -----------------------------------------------------------------

class ScoresPlotter(Plotter):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ScoresPlotter, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Get the data
        self.get_data()

        # 3. Plot
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

            x = []
            y = []
            ymax = []
            ymin = []
            yerr_max = []
            yerr_min = []
            ystddev = []
            diff_y = []

            # Loop over the generations
            for generation in get_generations(self.database, run_id):

                # Get the statistics
                statistics = get_statistics(self.database, run_id, generation)

                # GENERATION INDEX
                x.append(generation)

                # AVERAGE SCORE
                if self.config.fitness: y.append(statistics.fitness.average)
                else: y.append(statistics.raw.average)

                # MIN AND MAX
                if self.config.fitness:
                    ymax.append(statistics.fitness.max)
                    ymin.append(statistics.fitness.min)
                else:
                    ymax.append(statistics.raw.max)
                    ymin.append(statistics.raw.min)

                # SCORE ERROR BAR
                if self.config.fitness:
                    y_max = statistics.fitness.max - statistics.fitness.average
                    y_min = statistics.fitness.average - statistics.fitness.min
                else:
                    y_max = statistics.raw.max - statistics.raw.average
                    y_min = statistics.raw.average - statistics.raw.min
                yerr_max.append(y_max)
                yerr_min.append(y_min)

                # SCORE DIFFERENCE
                if self.config.fitness: diff_y.append(statistics.fitness.max - statistics.fitness.min)
                else: diff_y.append(statistics.raw.max - statistics.raw.min)

                # SCORE STDDEV
                if not self.config.fitness: ystddev.append(statistics.raw.stddev)

            if self.config.fitness: ystddev = None

            data = Map()
            data.x = x
            data.y = y
            data.ymax = ymax
            data.ymin = ymin
            data.yerr_max = yerr_max
            data.yerr_min = yerr_min
            data.ystddev = ystddev
            data.diff_y = diff_y

            # set the data
            self.data[run_id] = data

    # -----------------------------------------------------------------

    def plot(self):

        """
        THis function ...
        :return: 
        """

        # Inform the user
        log.info("Plotting ...")

        # Loop over the runs
        for run_id in self.runs: self.plot_scores(run_id)

        # Plot differences
        for run_id in self.runs: self.plot_differences(run_id)

        # Plot min/max/stddev
        for run_id in self.runs: self.plot_min_max(run_id)

    # -----------------------------------------------------------------

    def plot_scores(self, run_id):

        """
        This function ...
        :param run_id:
        :return: 
        """

        # Inform the user
        log.info("Plotting the scores for run '" + run_id + "' ...")

        plt.figure()

        x = self.data[run_id].x
        y = self.data[run_id].y
        yerr_min = self.data[run_id].yerr_min
        yerr_max = self.data[run_id].yerr_max

        plt.errorbar(x, y, yerr_min, yerr_max, ecolor="g")

        plt.xlabel("Generation")
        plt.ylabel("Score Min/Avg/Max")

        # Set title
        title = "Plot of evolution of run '" + run_id + "'"
        self.set_title(title)

        plt.grid(True)

        #plt.show()

        self.write_or_show("scores")

    # -----------------------------------------------------------------

    def plot_differences(self, run_id):

        """
        This function ...
        :param run_id:
        :return: 
        """

        # Inform the user
        log.info("Plotting the differences ...")

        plt.figure()

        x = self.data[run_id].x
        diff_y = self.data[run_id].diff_y

        plt.plot(x, diff_y, "g", label="Raw difference", linewidth=1.2)
        plt.fill_between(x, diff_y, color="g", alpha=0.1)

        diff_raw_max = max(diff_y)
        gen_max_raw = x[diff_y.index(diff_raw_max)]

        plt.annotate("Maximum (%.2f)" % (diff_raw_max,), xy=(gen_max_raw, diff_raw_max), xycoords='data',
                       xytext=(-150, -20), textcoords='offset points',
                       arrowprops=dict(arrowstyle="->",
                                       connectionstyle="arc"),
                       )

        plt.xlabel("Generation")
        plt.ylabel("Raw difference")
        self.set_title("Plot of score differences across evolution for run '" + run_id + "'")

        plt.grid(True)
        #plt.legend(prop=FontProperties(size="smaller"))

        #plt.show()
        self.write_or_show("differences")

    # -----------------------------------------------------------------

    def plot_min_max(self, run_id):

        """
        This function ...
        :param run_id: 
        :return: 
        """

        # Inform the user
        log.info("Plotting the minimal and maximal scores ...")

        plt.figure()

        x = self.data[run_id].x

        avg_y = self.data[run_id].y

        max_y = self.data[run_id].ymax
        min_y = self.data[run_id].ymin

        std_dev_y = self.data[run_id].ystddev

        plt.plot(x, max_y, "g", label="Max", linewidth=1.2)
        plt.plot(x, min_y, "r", label="Min", linewidth=1.2)
        plt.plot(x, avg_y, "b", label="Avg", linewidth=1.2)

        if std_dev_y is not None: plt.plot(x, std_dev_y, "k", label="Stddev", linewidth=1.2)

        plt.fill_between(x, min_y, max_y, color="g", alpha=0.1, label="Diff max/min")

        if self.config.minmax == "min": raw_max = min(min_y)
        else: raw_max = max(max_y)

        if self.config.minmax == "min": gen_max = x[min_y.index(raw_max)]
        else: gen_max = x[max_y.index(raw_max)]

        if std_dev_y is not None:
            min_std = min(std_dev_y)
            gen_min_std = x[std_dev_y.index(min_std)]

            max_std = max(std_dev_y)
            gen_max_std = x[std_dev_y.index(max_std)]

        if self.config.minmax == "min": annot_label = "Minimum (%.2f)" % (raw_max,)
        else: annot_label = "Maximum (%.2f)" % (raw_max,)

        plt.annotate(annot_label, xy=(gen_max, raw_max), xycoords='data',
                       xytext=(8, 15), textcoords='offset points',
                       arrowprops=dict(arrowstyle="->",
                                       connectionstyle="arc"),
                       )

        if std_dev_y is not None:

            plt.annotate("Min StdDev (%.2f)" % (min_std,), xy=(gen_min_std, min_std), xycoords='data',
                           xytext=(8, 15), textcoords='offset points',
                           arrowprops=dict(arrowstyle="->",
                                           connectionstyle="arc"),
                           )

            plt.annotate("Max StdDev (%.2f)" % (max_std,), xy=(gen_max_std, max_std), xycoords='data',
                           xytext=(8, 15), textcoords='offset points',
                           arrowprops=dict(arrowstyle="->",
                                           connectionstyle="arc"),
                           )

        plt.xlabel("Generation")
        plt.ylabel("Raw score")

        self.set_title("Plot of evolution for run '" + run_id + "'")

        plt.grid(True)
        #plt.legend(prop=FontProperties(size="smaller"))

        #plt.show()

        self.write_or_show("minmax")

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
