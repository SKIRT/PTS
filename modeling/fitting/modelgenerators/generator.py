#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.modelgenerators.generator Contains the abstract ModelGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....core.tools.logging import log
from ..component import FittingComponent

# -----------------------------------------------------------------

class ModelGenerator(FittingComponent):
    
    """
    This class...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(ModelGenerator, self).__init__()

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the necessary input
        self.load_input()

        # Initialize the animations
        self.initialize_animations()

        # 2. Generate the model parameters
        self.generate()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup of the base class
        super(ModelGenerator, self).setup()

    # -----------------------------------------------------------------

    def load_input(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def initialize_animations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Initializing the animations ...")

        # Initialize the scatter animation
        self.scatter_animation = ScatterAnimation(self.ranges["FUV young"], self.ranges["FUV ionizing"],
                                                  self.ranges["Dust mass"])
        self.scatter_animation.x_label = "FUV luminosity of young stars"
        self.scatter_animation.y_label = "FUV luminosity of ionizing stars"
        self.scatter_animation.z_label = "Dust mass"

        # Initialize the young FUV luminosity distribution animation
        self.fuv_young_animation = DistributionAnimation(self.ranges["FUV young"][0], self.ranges["FUV young"][1],
                                                         "FUV luminosity of young stars", "New models")
        self.fuv_young_animation.add_reference_distribution("Previous models", self.distributions["FUV young"])

        # Initialize the ionizing FUV luminosity distribution animation
        self.fuv_ionizing_animation = DistributionAnimation(self.ranges["FUV ionizing"][0],
                                                            self.ranges["FUV ionizing"][1],
                                                            "FUV luminosity of ionizing stars", "New models")
        self.fuv_ionizing_animation.add_reference_distribution("Previous models", self.distributions["FUV ionizing"])

        # Initialize the dust mass distribution animation
        self.dust_mass_animation = DistributionAnimation(self.ranges["Dust mass"][0], self.ranges["Dust mass"][1],
                                                         "Dust mass", "New models")
        self.dust_mass_animation.add_reference_distribution("Previous models", self.distributions["Dust mass"])

    # -----------------------------------------------------------------

    def generate(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def update_animations(self, young_luminosity, ionizing_luminosity, dust_mass):

        """
        This function ...
        :param young_luminosity:
        :param ionizing_luminosity:
        :param dust_mass:
        :return:
        """

        # Add the point (and thus a frame) to the animation of parameter points
        self.scatter_animation.add_point(young_luminosity, ionizing_luminosity, dust_mass)

        # Update the distribution animations
        if self.number_of_models > 1:

            # Add a frame to the animation of the distribution of the FUV luminosity of young starss
            self.fuv_young_animation.add_value(young_luminosity)

            # Add a frame to the animation of the distribution of the FUV luminosity of ionizing stars
            self.fuv_ionizing_animation.add_value(ionizing_luminosity)

            # Add a frame to the animation of the distribution of the dust mass
            self.dust_mass_animation.add_value(dust_mass)

    # -----------------------------------------------------------------

    @property
    def number_of_models(self):

        """
        This function ...
        :return:
        """

        return len(self.parameters["FUV young"])

# -----------------------------------------------------------------
