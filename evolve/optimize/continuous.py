#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.evolve.optimize.continuous Contains the ContinuousOptimizer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from .optimizer import Optimizer, show_best

# -----------------------------------------------------------------

class ContinuousOptimizer(Optimizer):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ContinuousOptimizer, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Initialize
        self.initialize(**kwargs)

        # 3. Evolve
        self.evolve()

        # 4. Show
        if self.config.show: self.show()

        # 5. Write
        if self.config.write: self.write()

        # 6. Plot
        if self.config.plot: self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ContinuousOptimizer, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def initialize(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Inform the user
        log.info("Initializing ...")

        # 1. Initialize the random number generator
        self.initialize_prng()

        # 2. Initialize the adapters
        self.initialize_adapters()

        # 3. Initialize genome
        self.initialize_genome(**kwargs)

        # 4. Set the evaluator of the initial genome
        self.initial_genome.evaluator.set(kwargs.pop("evaluator"))

        # 5. Initialize engine
        self.initialize_engine(**kwargs)

        # 6. Initialize population
        self.initialize_population()

    # -----------------------------------------------------------------

    def evolve(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Evolving ...")

        # Set the database adapters again # NO: NOT NECESSARY IN CONTINUOUS OPTIMIZER: INITIALIZE() HAS DONE THIS IN INITIALIZE_ENGINE, AND THIS CLASS STAYS IN CURRENT RUNTIME
        #self.set_engine_adapters()

        # Let evolve
        self.engine.evolve(freq_stats=self.config.stats_freq, progress_bar=(not log.is_debug))

        # Get the best individual
        self.best = self.engine.bestIndividual()

        # Plot best
        if self.generations_plotter is not None: self.generations_plotter.add_best(self.best)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

        # Show the best individual
        show_best(self.best)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # 1. Write the database
        self.write_database()

        # 2. Write the statistics
        self.write_statistics()

        # 3. Write the best individual
        self.write_best()

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot scores
        self.plot_scores()

        # Plot heat map
        self.plot_heat_map()

    # -----------------------------------------------------------------

    def plot_scores(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the scores ...")

    # -----------------------------------------------------------------

    def plot_heat_map(self):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------
