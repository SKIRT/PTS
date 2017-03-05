#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.modelgenerators.instinctive Contains the InstinctiveModelGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....core.tools.logging import log
from .generator import ModelGenerator

# -----------------------------------------------------------------

class InstinctiveModelGenerator(ModelGenerator):
    
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
        super(InstinctiveModelGenerator, self).__init__(config, interactive)

        # The fitting run
        self.fitting_run = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(InstinctiveModelGenerator, self).setup(**kwargs)

        # Get the fitting run
        self.fitting_run = kwargs.pop("fitting_run")

    # -----------------------------------------------------------------

    def load_input(self):

        """
        This function ...
        :return:
        """

        self.load_distributions()

    # -----------------------------------------------------------------

    def load_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the probability distributions for the different fit parameters ...")

        # Loop over the different fit parameters
        for parameter_name in self.parameter_names:

            # Load the probability distribution
            distribution = Distribution.from_file(self.distribution_table_paths[parameter_name])

            # Normalize the distribution
            distribution.normalize(value=1.0, method="max")

            # Set the distribution
            self.distributions[parameter_name] = distribution

    # -----------------------------------------------------------------

    def generate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Picking random parameter values based on the probability distributions ...")

        # Draw parameters values for the specified number of simulations
        for counter in range(self.config.simulations):

            # Debugging
            log.debug("Calculating random parameter set " + str(counter + 1) + " of " + str(self.config.simulations) + " ...")

            # Draw a random FUV luminosity of the young stellar population
            young_luminosity = self.distributions["FUV young"].random(self.config.young_stars.min, self.config.young_stars.max)

            # Draw a random FUV luminosity of the ionizing stellar population
            ionizing_luminosity = self.distributions["FUV ionizing"].random(self.config.ionizing_stars.min, self.config.ionizing_stars.max)

            # Draw a random dust mass
            dust_mass = self.distributions["Dust mass"].random(self.config.dust.min, self.config.dust.max)

            # Add the parameter values to the dictionary
            self.parameters["FUV young"].append(young_luminosity)
            self.parameters["FUV ionizing"].append(ionizing_luminosity)
            self.parameters["Dust mass"].append(dust_mass)

            # Update the animations
            if self.config.visualise: self.update_animations(young_luminosity, ionizing_luminosity, dust_mass)

# -----------------------------------------------------------------
