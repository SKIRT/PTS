#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.modelgenerators.grid Contains the GridModelGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....core.tools.logging import log
from .generator import ModelGenerator

# -----------------------------------------------------------------

class GridModelGenerator(ModelGenerator):
    
    """
    This class...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(GridModelGenerator, self).__init__(config, interactive)

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the constructor of the base class
        super(GridModelGenerator, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def generate(self):

        """
        THis function ...
        :return: 
        """

        # Inform the user
        log.info("")

        scale = None
        grid_points = self.generate_grid_points(scale) # returns dictionary

        # FROM GENETIC:

        # NAMED:

        # Loop over the individual names
        for name in self.individual_names:

            # Get the individual
            individual = self.optimizer.population[name]

            # Loop over all the genes (parameters)
            for i in range(len(individual)):

                # Get the parameter value
                value = individual[i]

                # Add the parameter value to the dictionary
                self.parameters[self.fitting_run.free_parameter_labels[i]][name] = value

        # UNNAMED:

        # Loop over the individuals of the population
        for individual in self.optimizer.population:

            # Loop over all the genes (parameters)
            for i in range(len(individual)):
                # Get the parameter value
                value = individual[i]

                # Add the parameter value to the dictionary
                self.parameters[self.fitting_run.free_parameter_labels[i]].append(value)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return: 
        """

# -----------------------------------------------------------------
