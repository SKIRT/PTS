#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.evolve.solve.extremizer Contains the Extremizer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable

# -----------------------------------------------------------------

def find_minimum(function, initial_values, ranges):

    """
    This function ...
    :param function:
    :param initial_values:
    :param ranges:
    :return: 
    """

# -----------------------------------------------------------------

def find_maximum(function, initial_values, ranges):

    """
    This function ...
    :param function:
    :return: 
    """

# -----------------------------------------------------------------

class Extremizer(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None, interactive=False):

        """
        This function ...
        :param config: 
        :param interactive: 
        """

        # Call the constructor of the base class
        super(Extremizer, self).__init__(config, interactive)

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs: 
        :return: 
        """

        # Call the setup function
        self.setup(**kwargs)

        # Optimize
        self.optimize()

        # Plot
        self.plot()

        # Write
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs: 
        :return: 
        """

        # Call the setup function of the base class
        super(Extremizer, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def optimize(self):

        """
        This function ...
        :return: 
        """

        # Settings
        settings_optimize = dict()
        settings_optimize["output"] = None
        settings_optimize["nparameters"] = nparameters
        settings_optimize["nindividuals"] = nindividuals
        settings_optimize["parameter_range"] = parameter_range
        settings_optimize["best_raw_score"] = best_raw_score
        settings_optimize["round_decimal"] = round_decimal
        settings_optimize["ngenerations"] = ngenerations
        settings_optimize["mutation_rate"] = mutation_rate
        settings_optimize["crossover_rate"] = crossover_rate
        settings_optimize["stats_freq"] = stats_freq
        settings_optimize["mutation_method"] = mutation_method
        settings_optimize["min_or_max"] = min_or_max

        # Input
        input_optimize = dict()
        input_optimize["evaluator"] = rastringin

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Plotting ...")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return: 
        """

# -----------------------------------------------------------------
