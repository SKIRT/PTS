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

# Import standard modules
from abc import abstractmethod, ABCMeta
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ....core.tools.logging import log
from ..component import FittingComponent
from ....magic.animation.scatter import ScatterAnimation
from ....magic.animation.distribution import DistributionAnimation

# -----------------------------------------------------------------

class ModelGenerator(FittingComponent):
    
    """
    This class...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(ModelGenerator, self).__init__()

        # The dictionary with the parameter ranges
        self.ranges = OrderedDict()

        # The dictionary with the list of the model parameters
        self.parameters = OrderedDict()

        # The animations
        self.scatter_animation = None
        self.parameter_animations = dict()

    # -----------------------------------------------------------------

    @property
    def parameter_names(self):

        """
        This function ...
        :return:
        """

        return self.ranges.keys()

    # -----------------------------------------------------------------

    @property
    def nparameters(self):

        """
        This function ...
        :return:
        """

        return len(self.ranges)

    # -----------------------------------------------------------------

    @property
    def nmodels(self):

        """
        This function ...
        :return:
        """

        return len(self.parameters[self.ranges.keys()[0]])

    # -----------------------------------------------------------------

    @property
    def parameter_minima(self):

        """
        This function ...
        :return:
        """

        minima = []
        for name in self.ranges: minima.append(self.ranges[name].min)

        # Return the minimal parameter values
        return minima

    # -----------------------------------------------------------------

    @property
    def parameter_maxima(self):

        """
        This function ...
        :return:
        """

        maxima = []
        for name in self.ranges: maxima.append(self.ranges[name].max)

        # Return the maximal parameter values
        return maxima

    # -----------------------------------------------------------------

    def add_parameter(self, name, par_range):

        """
        This function ...
        :param name:
        :param par_range:
        :return:
        """

        self.ranges[name] = par_range

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

        # 3. Initialize the animations
        self.initialize_animations()

        # 4. Generate the model parameters
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

        # Initialize the scatter animation, if there are exactly 3 free parameters
        if len(self.free_parameter_labels) == 3:

            self.scatter_animation = ScatterAnimation(self.ranges["FUV young"], self.ranges["FUV ionizing"], self.ranges["Dust mass"])
            self.scatter_animation.x_label = "FUV luminosity of young stars"
            self.scatter_animation.y_label = "FUV luminosity of ionizing stars"
            self.scatter_animation.z_label = "Dust mass"

        # Loop over the free parameters and create an individual animation for each of them
        for label in self.free_parameter_labels:

            description = self.free_parameter_labels[label]

            # Initialize the young FUV luminosity distribution animation
            animation = DistributionAnimation(self.ranges[label][0], self.ranges[label][1], description, "New models")
            animation.add_reference_distribution("Previous models", self.distributions[label])

            # Add the animation to the dictionary
            self.parameter_animations[label] = animation

    # -----------------------------------------------------------------

    @abstractmethod
    def generate(self):

        """
        This function ...
        :return:
        """

        pass

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
        if self.nmodels > 1:

            # Add a frame to the animation of the distribution of the FUV luminosity of young starss
            self.fuv_young_animation.add_value(young_luminosity)

            # Add a frame to the animation of the distribution of the FUV luminosity of ionizing stars
            self.fuv_ionizing_animation.add_value(ionizing_luminosity)

            # Add a frame to the animation of the distribution of the dust mass
            self.dust_mass_animation.add_value(dust_mass)

# -----------------------------------------------------------------
