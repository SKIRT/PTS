#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.geneticevolver Contains the GeneticModelEvolver class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .evolver import ModelEvolver
from ...core.tools.logging import log
from ...magic.animation.scatter import ScatterAnimation
from ...magic.animation.distribution import DistributionAnimation

# -----------------------------------------------------------------

class GeneticModelEvolver(ModelEvolver):
    
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
        super(GeneticModelEvolver, self).__init__(config)

        # -- Attributes --

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new ParameterExplorer instance
        explorer = cls(arguments.config)

        # Set the modeling path
        explorer.config.path = arguments.path

        # Set the number of simulations to launch in the batch
        if arguments.simulations is not None: explorer.config.simulations = arguments.simulations

        # Set the remote host IDs
        if arguments.remotes is not None: explorer.config.remotes = arguments.remotes

        # Set the limits of the FUV luminosity of the young stellar population
        if arguments.young is not None:
            explorer.config.young_stars.min = arguments.young[0]
            explorer.config.young_stars_max = arguments.young[1]

        # Set the limits of the FUV luminosity of the ionizing stellar population
        if arguments.ionizing is not None:
            explorer.config.ionizing_stars.min = arguments.ionizing[0]
            explorer.config.ionizing_stars.max = arguments.ionizing[1]

        # Set the limits of the dust mass
        if arguments.dust is not None:
            explorer.config.dust.min = arguments.dust[0]
            explorer.config.dust.max = arguments.dust[1]

        # Make visualisations
        explorer.config.visualise = arguments.visualise

        # Return the new instance
        return explorer

    # -----------------------------------------------------------------

    def initialize_animations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Initializing the animations ...")

        # Initialize the scatter animation
        self.scatter_animation = ScatterAnimation([self.config.young_stars.min, self.config.young_stars.max],
                                                  [self.config.ionizing_stars.min, self.config.ionizing_stars.max],
                                                  [self.config.dust.min, self.config.dust.max])
        self.scatter_animation.x_label = "FUV luminosity of young stars"
        self.scatter_animation.y_label = "FUV luminosity of ionizing stars"
        self.scatter_animation.z_label = "Dust mass"

        # Initialize the young FUV luminosity distribution animation
        self.fuv_young_animation = DistributionAnimation(self.config.young_stars.min, self.config.young_stars.max, "FUV luminosity of young stars", "New models")

        # Initialize the ionizing FUV luminosity distribution animation
        self.fuv_ionizing_animation = DistributionAnimation(self.config.ionizing_stars.min, self.config.ionizing_stars.max, "FUV luminosity of ionizing stars", "New models")

        # Initialize the dust mass distribution animation
        self.dust_mass_animation = DistributionAnimation(self.config.dust.min, self.config.dust.max, "Dust mass", "New models")

    # -----------------------------------------------------------------

    def set_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Picking new individuals for the next generations ...")

        # Draw parameters values for the specified number of simulations
        #for counter in range(self.config.simulations):


# -----------------------------------------------------------------
