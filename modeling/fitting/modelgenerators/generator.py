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
from collections import OrderedDict, defaultdict

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

        # The order of the free parameter labels
        self.parameter_labels_order = None

        # The dictionary with the parameter ranges
        self.ranges = OrderedDict()

        # The dictionary with the list of the model parameters
        self.parameters = defaultdict(list)

        # The parameter value distributions
        self.distributions = dict()

        # The animations
        self.scatter_animation = None
        self.scatter_animation_labels = None
        self.parameter_animations = dict()

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

        return len(self.parameters[self.parameter_labels_order[0]])

    # -----------------------------------------------------------------

    @property
    def parameter_minima(self):

        """
        This function ...
        :return:
        """

        # Initialize a list
        minima = []

        # Set the list values
        for label in self.parameter_labels_order: minima.append(self.ranges[label].min)

        # Return the minimal parameter values
        return minima

    # -----------------------------------------------------------------

    @property
    def parameter_maxima(self):

        """
        This function ...
        :return:
        """

        # Initialize a list
        maxima = []

        # Set the list values
        for label in self.parameter_labels_order: maxima.append(self.ranges[label].max)

        # Return the maximal parameter values
        return maxima

    # -----------------------------------------------------------------

    def add_parameter(self, label, parameter_range):

        """
        This function ...
        :param label:
        :param parameter_range:
        :return:
        """

        self.ranges[label] = parameter_range

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 3. Load the current parameter value probability distributions
        self.load_distributions()

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

        # Establish the order of the free parameters
        self.parameter_labels_order = list(self.free_parameter_labels)

    # -----------------------------------------------------------------

    def load_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the current parameter distributions ...")

        # Loop over the free parameters
        for label in self.free_parameter_labels:

            # Load the distribution
            if self.has_distribution(label): self.distributions[label] = self.get_parameter_distribution(label)

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

            label0, label1, label2 = self.free_parameter_labels
            description0, description1, description2 = [self.parameter_descriptions[label] for label in [label0, label1, label2]]

            # Establish the order of the free parameters in the scatter ani
            self.scatter_animation_labels = [label0, label1, label2]

            # Initialize the scatter animation
            self.scatter_animation = ScatterAnimation(self.ranges[label0], self.ranges[label1], self.ranges[label2])
            self.scatter_animation.x_label = description0
            self.scatter_animation.y_label = description1
            self.scatter_animation.z_label = description2

        # Loop over the free parameters and create an individual animation for each of them
        for label in self.free_parameter_labels:

            description = self.parameter_descriptions[label]

            # Initialize the young FUV luminosity distribution animation
            animation = DistributionAnimation(self.ranges[label].min, self.ranges[label].max, description, "New models")
            if label in self.distributions: animation.add_reference_distribution("Previous models", self.distributions[label])

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

    def update_animations(self, values_dict):

        """
        This function ...
        :param values_dict:
        :return:
        """

        # Add the point (and thus a frame) to the animation of parameter points
        if self.scatter_animation is not None:

            # Set the values
            values = []
            for label in self.scatter_animation_labels:
                values.append(values_dict[label])

            self.scatter_animation.add_point(values[0], values[1], values[2])

        # Update the distribution animations
        if self.nmodels > 1:

            # Add a point (and thus one frame) to the individual parameter animations
            for label in values_dict: self.parameter_animations[label].add_value(values_dict[label])

# -----------------------------------------------------------------
