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
from ...core.tools.logging import log
from .optimizer import Optimizer, show_best

# -----------------------------------------------------------------

class ContinuousOptimizer(Optimizer):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        This function ...
        :param config:
        """

        # Call the constructor of the base class
        super(ContinuousOptimizer, self).__init__(config)

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

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

        # 1. Initialize the database
        self.initialize_database()

        # 2. Initialize genome
        self.initialize_genome(**kwargs)

        # 3. Set the evaluator of the initial genome
        self.initial_genome.evaluator.set(kwargs.pop("evaluator"))

        # 4. Initialize engine
        self.initialize_engine(**kwargs)

    # -----------------------------------------------------------------

    def evolve(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Evolving ...")

        # Let evolve
        self.engine.evolve(freq_stats=self.config.stats_freq, progress_bar=self.config.progress_bar)

        # Get the best individual
        self.best = self.engine.bestIndividual()

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

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

# -----------------------------------------------------------------
