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
from ...core.basics.configurable import Configurable
from ...core.tools import filesystem as fs
from .plotter import Plotter
from ...core.tools.logging import log
from ..analyse.database import get_individuals, get_generations, get_runs

# -----------------------------------------------------------------

class ScoresPlotter(Plotter):

    """
    This class ...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        """

        # Call the constructor of the base class
        super(ScoresPlotter, self).__init__(config, interactive)

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # Plot
        self.plot_scores()

    # -----------------------------------------------------------------

    @property
    def runs(self):

        """
        This function ...
        :return: 
        """

        if self.config.runs is None: get_runs(self.database)
        else: return self.config.runs

    # -----------------------------------------------------------------

    def plot_scores(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Plotting ...")

        x = []
        y = []
        yerr_max = []
        yerr_min = []

        # Loop over the runs
        for run_id in self.runs:

            # Loop over the generations
            for generation in get_generations(self.database, run_id):

                # Loop over the individuals
                for ind in get_individuals(self.database, run_id, generation):

                    # Get the name
                    #name = ind[""]

                    x.append(ind["generation"])
                    y.append(ind["rawAve"])
                    ymax = ind["rawMax"] - ind["rawAve"]
                    ymin = ind["rawAve"] - ind["rawMin"]

                    yerr_max.append(ymax)
                    yerr_min.append(ymin)

        plt.figure()

        plt.errorbar(x, y, [yerr_min, yerr_max], ecolor="g")

        plt.xlabel('Generation (#)')
        plt.ylabel('Raw score Min/Avg/Max')
        #plt.title("Plot of evolution identified by '%s' (raw scores)" % (options.identify))
        plt.grid(True)

        plt.show()

    # -----------------------------------------------------------------

    def graph_errorbars_raw(self, pop, minimize, filesave=None):

        """
        This function ...
        :param pop:
        :param minimize:
        :param filesave:
        :return:
        """

        x = []
        y = []
        yerr_max = []
        yerr_min = []

        for it in pop:
            x.append(it["generation"])
            y.append(it["rawAve"])
            ymax = it["rawMax"] - it["rawAve"]
            ymin = it["rawAve"] - it["rawMin"]

            yerr_max.append(ymax)
            yerr_min.append(ymin)

        pylab.figure()
        pylab.errorbar(x, y, [yerr_min, yerr_max], ecolor="g")
        pylab.xlabel('Generation (#)')
        pylab.ylabel('Raw score Min/Avg/Max')
        pylab.title("Plot of evolution identified by '%s' (raw scores)" % (options.identify))
        pylab.grid(True)

        if filesave:
            pylab.savefig(filesave)
            print("Graph saved to %s file !" % (filesave,))
        else:
            pylab.show()

    # -----------------------------------------------------------------

    def graph_diff_raw(self, pop, minimize, filesave=None):

        """
        This function ...
        :param pop:
        :param minimize:
        :param filesave:
        :return:
        """

        x = []

        diff_raw_y = []
        diff_fit_y = []

        for it in pop:

            x.append(it["generation"])
            diff_raw_y.append(it["rawMax"] - it["rawMin"])
            diff_fit_y.append(it["fitMax"] - it["fitMin"])

        pylab.figure()
        pylab.subplot(211)

        pylab.plot(x, diff_raw_y, "g", label="Raw difference", linewidth=1.2)
        pylab.fill_between(x, diff_raw_y, color="g", alpha=0.1)

        diff_raw_max = max(diff_raw_y)
        gen_max_raw = x[diff_raw_y.index(diff_raw_max)]

        pylab.annotate("Maximum (%.2f)" % (diff_raw_max,), xy=(gen_max_raw, diff_raw_max), xycoords='data',
                       xytext=(-150, -20), textcoords='offset points',
                       arrowprops=dict(arrowstyle="->",
                                       connectionstyle="arc"),
                       )

        pylab.xlabel("Generation (#)")
        pylab.ylabel("Raw difference")
        pylab.title("Plot of evolution identified by '%s'" % (options.identify))

        pylab.grid(True)
        pylab.legend(prop=FontProperties(size="smaller"))

        pylab.subplot(212)

        pylab.plot(x, diff_fit_y, "b", label="Fitness difference", linewidth=1.2)
        pylab.fill_between(x, diff_fit_y, color="b", alpha=0.1)

        diff_fit_max = max(diff_fit_y)
        gen_max_fit = x[diff_fit_y.index(diff_fit_max)]

        pylab.annotate("Maximum (%.2f)" % (diff_fit_max,), xy=(gen_max_fit, diff_fit_max), xycoords='data',
                       xytext=(-150, -20), textcoords='offset points',
                       arrowprops=dict(arrowstyle="->",
                                       connectionstyle="arc"),)

        pylab.xlabel("Generation (#)")
        pylab.ylabel("Fitness difference")

        pylab.grid(True)
        pylab.legend(prop=FontProperties(size="smaller"))

        if filesave:
            pylab.savefig(filesave)
            print("Graph saved to %s file !" % (filesave,))
        else: pylab.show()

# -----------------------------------------------------------------
