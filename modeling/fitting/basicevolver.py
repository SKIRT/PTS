#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.advancedparameterexplorer Contains the AdvancedParameterExplorer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .evolver import ModelEvolver
from ...core.tools.logging import log
from ...core.basics.distribution import Distribution
from ...magic.animation.scatter import ScatterAnimation
from ...magic.animation.distribution import DistributionAnimation

# -----------------------------------------------------------------

class BasicModelEvolver(ModelEvolver):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(BasicModelEvolver, self).__init__(config)

        # -- Attributes --

        # The probability distributions for the different fit parameters
        self.distributions = dict()

    # -----------------------------------------------------------------

    def load_input(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the necessary input ...")

        # 1. Load the current parameter table
        self.load_table()

        # 2. Load the ski file
        self.load_ski()

        # 3. Load the probability distributions for the different parameters
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

    def initialize_animations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Initializing the animations ...")

        # Initialize the scatter animation
        self.scatter_animation = ScatterAnimation(self.ranges["FUV young"], self.ranges["FUV ionizing"], self.ranges["Dust mass"])
        self.scatter_animation.x_label = "FUV luminosity of young stars"
        self.scatter_animation.y_label = "FUV luminosity of ionizing stars"
        self.scatter_animation.z_label = "Dust mass"

        # Initialize the young FUV luminosity distribution animation
        self.fuv_young_animation = DistributionAnimation(self.ranges["FUV young"][0], self.ranges["FUV young"][1], "FUV luminosity of young stars", "New models")
        self.fuv_young_animation.add_reference_distribution("Previous models", self.distributions["FUV young"])

        # Initialize the ionizing FUV luminosity distribution animation
        self.fuv_ionizing_animation = DistributionAnimation(self.ranges["FUV ionizing"][0], self.ranges["FUV ionizing"][1], "FUV luminosity of ionizing stars", "New models")
        self.fuv_ionizing_animation.add_reference_distribution("Previous models", self.distributions["FUV ionizing"])

        # Initialize the dust mass distribution animation
        self.dust_mass_animation = DistributionAnimation(self.ranges["Dust mass"][0], self.ranges["Dust mass"][1], "Dust mass", "New models")
        self.dust_mass_animation.add_reference_distribution("Previous models", self.distributions["Dust mass"])

    # -----------------------------------------------------------------

    def set_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Picking random parameter values based on the probability distributions ...")

        # Draw parameters values for the specified number of simulations
        for counter in range(self.config.simulations):

            # Debugging
            log.debug("Calculating random parameter set " + str(counter+1) + " of " + str(self.config.simulations) + " ...")

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
