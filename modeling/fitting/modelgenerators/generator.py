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
from abc import abstractmethod, ABCMeta, abstractproperty
from collections import OrderedDict, defaultdict

# Import the relevant PTS classes and modules
from ....core.basics.log import log
from ..component import FittingComponent
from ....magic.animation.scatter import ScatterAnimation
from ....magic.animation.distribution import DistributionAnimation
from ....core.tools import types
from ....core.tools import numbers, sequences

# -----------------------------------------------------------------

class ModelGenerator(FittingComponent):
    
    """
    This class...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(ModelGenerator, self).__init__(*args, **kwargs)

        # The fitting run
        self.fitting_run = None

        # The generation
        self.generation = None

        # The dictionary with the parameter ranges
        self.ranges = OrderedDict()

        # The dictionary with the list of the model parameters
        self.parameters = defaultdict(dict)

        # The parameter value distributions
        self.distributions = dict()

        # The animations
        self.scatter_animation = None
        self.scatter_animation_labels = None
        self.parameter_animations = dict()

        # The scales for the different parameters
        self.scales = None

        # Other input
        self.most_sampled_parameters = None
        self.sampling_weights = None

    # -----------------------------------------------------------------

    @property
    def generation_path(self):

        """
        This function ...
        :return:
        """

        return self.generation.path

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

        return len(self.parameters[self.fitting_run.free_parameter_labels[0]])

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
        for label in self.fitting_run.free_parameter_labels: minima.append(self.ranges[label].min)

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
        for label in self.fitting_run.free_parameter_labels: maxima.append(self.ranges[label].max)

        # Return the maximal parameter values
        return maxima

    # -----------------------------------------------------------------

    @property
    def parameter_minima_scalar(self):

        """
        This function ...
        :return:
        """

        # Initialize a list
        minima = []

        # Set the list values
        for label in self.fitting_run.free_parameter_labels:

            min_value = self.ranges[label].min

            # Convert if necessary
            if label in self.fitting_run.parameter_units and self.fitting_run.parameter_units[label] is not None:
                unit = self.fitting_run.parameter_units[label]
                min_value = min_value.to(unit).value

            # Assert that is real type
            assert types.is_real_type(min_value)
            min_value = float(min_value)

            # Add to list
            minima.append(min_value)

        # Return the minimal parameter values
        return minima

    # -----------------------------------------------------------------

    @property
    def parameter_maxima_scalar(self):

        """
        This function ...
        :return:
        """

        # Initialize a list
        maxima = []

        # Set the list values
        for label in self.fitting_run.free_parameter_labels:

            max_value = self.ranges[label].max

            # Convert if necessary
            if label in self.fitting_run.parameter_units and self.fitting_run.parameter_units[label] is not None:
                unit = self.fitting_run.parameter_units[label]
                max_value = max_value.to(unit).value

            # Assert that is real type
            assert types.is_real_type(max_value)
            max_value = float(max_value)

            # Add to list
            maxima.append(max_value)

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

        # Add range
        self.ranges[label] = parameter_range

    # -----------------------------------------------------------------

    @abstractproperty
    def individual_names(self):

        """
        This function ...
        :return: 
        """

        pass

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 3. Load the current parameter value probability distributions
        self.load_distributions()

        # 3. Initialize the animations
        if self.config.animate: self.initialize_animations()

        # 4. Generate the model parameters
        self.generate()

        # 5. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup of the base class
        super(ModelGenerator, self).setup(**kwargs)

        # Get the fitting run
        self.fitting_run = kwargs.pop("fitting_run")

        # Get the generation
        self.generation = kwargs.pop("generation")

        # Get scales for different free parameters
        if "scales" in kwargs: self.scales = kwargs.pop("scales")
        else: self.scales = self.config.scales
        #if self.scales is None: raise ValueError("Scales are not defined")

        # Get parameter ranges
        if "parameter_ranges" in kwargs: self.ranges = kwargs.pop("parameter_ranges")

        # Get other input
        if "most_sampled_parameters" in kwargs: self.most_sampled_parameters = kwargs.pop("most_sampled_parameters")
        if "sampling_weights" in kwargs: self.sampling_weights = kwargs.pop("sampling_weights")

    # -----------------------------------------------------------------

    def load_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the current parameter distributions ...")

        # Loop over the free parameters
        for label in self.fitting_run.free_parameter_labels:

            # Load the distribution
            if self.fitting_run.has_distribution(label): self.distributions[label] = self.fitting_run.get_parameter_distribution(label)

    # -----------------------------------------------------------------

    def initialize_animations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Initializing the animations ...")

        # Initialize the scatter animation, if there are exactly 3 free parameters
        if len(self.fitting_run.free_parameter_labels) == 3:

            label0, label1, label2 = self.fitting_run.free_parameter_labels
            description0, description1, description2 = [self.fitting_run.parameter_descriptions[label] for label in [label0, label1, label2]]

            # Establish the order of the free parameters in the scatter ani
            self.scatter_animation_labels = [label0, label1, label2]

            # Initialize the scatter animation
            self.scatter_animation = ScatterAnimation(self.ranges[label0], self.ranges[label1], self.ranges[label2])
            self.scatter_animation.x_label = description0
            self.scatter_animation.y_label = description1
            self.scatter_animation.z_label = description2

        # Loop over the free parameters and create an individual animation for each of them
        for label in self.fitting_run.free_parameter_labels:

            description = self.fitting_run.parameter_descriptions[label]

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

    def generate_grid_points_one_scale(self, scale, npoints=None, most_sampled=None, weights=None):

        """
        This function ...
        :param scale:
        :param npoints:
        :param most_sampled:
        :param weights:
        :return: 
        """

        # Inform the user
        log.info("Generating grid points for the parameter values of models on a " + scale + " scale ...")

        # The order is important!
        if most_sampled is not None: sampled_most = [(label in most_sampled) for label in self.fitting_run.free_parameter_labels]
        else: sampled_most = None

        # If the number of points is not defined, determine based on the total number of models
        if npoints is None:

            # Initialize dictionary of the number of points per parameter
            npoints = dict()

            # Determine the number of grid points for each parameter based on the desired total number of models
            npoints_list = numbers.divide_in_n_dimensions(self.config.nmodels, self.fitting_run.nfree_parameters, sampled_most=sampled_most, weights=weights)
            for label_index, label in enumerate(self.fitting_run.free_parameter_labels): npoints[label] = npoints_list[label_index]

        # Initialize a dictionary for the grid points for each parameter
        grid_points = dict()

        # Generate the grid points for each parameter, and loop from the center of the range
        for label in self.fitting_run.free_parameter_labels:

            # Determine the number of grid points
            npoints_parameter = npoints[label]

            # Get the range
            parameter_range = self.ranges[label]

            # Generate the grid points, based on the scale
            values = generate_grid_points_from_center(parameter_range, npoints_parameter, scale)

            # Set the values
            grid_points[label] = values

        # Return the dictionary
        return grid_points

    # -----------------------------------------------------------------

    def generate_grid_points_different_scales(self, scales, npoints=None, most_sampled=None, weights=None):

        """
        This function ...
        :param scales:
        :param npoints:
        :param most_sampled:
        :param weights:
        :return: 
        """

        # Inform the user
        log.info("Generating grid points for the parameter values of models ...")

        # The order is important!
        if most_sampled is not None: sampled_most = [(label in most_sampled) for label in self.fitting_run.free_parameter_labels]
        else: sampled_most = None

        # If the number of points is not defined, determine based on the total number of models
        if npoints is None:

            # Initialize dictionary of the number of points per parameter
            npoints = dict()

            # Determine the number of grid points for each parameter based on the desired total number of models
            npoints_list = numbers.divide_in_n_dimensions(self.config.nmodels, self.fitting_run.nfree_parameters, sampled_most=sampled_most, weights=weights)
            for label_index, label in enumerate(self.fitting_run.free_parameter_labels): npoints[label] = npoints_list[label_index]

        # Initialize a dictionary for the grid points for each parameter
        grid_points = dict()

        # Generate the grid points for each parameter, and loop from the center of the range
        for label in self.fitting_run.free_parameter_labels:

            # Determine the number of grid points
            npoints_parameter = npoints[label]

            # Get the range
            parameter_range = self.ranges[label]

            # Get the scale
            scale = scales[label]

            # Generate the grid points, based on the scale
            values = generate_grid_points_from_center(parameter_range, npoints_parameter, scale)

            # Set the values
            grid_points[label] = values

        # Return the dictionary
        return grid_points

    # -----------------------------------------------------------------

    def grid_points_to_lists(self, grid_points_dict):

        """
        This function ...
        :param grid_points_dict:
        :return: 
        """

        # Convert into lists, and strip units
        grid_points_lists = []
        for label in self.fitting_run.free_parameter_labels:

            # Get the list of scalar values
            if label in self.fitting_run.parameter_units and self.fitting_run.parameter_units[label] is not None:
                unit = self.fitting_run.parameter_units[label]
                #print(label, unit)
                values = [value.to(unit).value for value in grid_points_dict[label]]
            else: values = grid_points_dict[label]

            # Add the list of grid point values
            grid_points_lists.append(values)

        # Return the lists
        return grid_points_lists

    # -----------------------------------------------------------------

    @abstractmethod
    def write(self):

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

def generate_grid_points_from_center(parameter_range, npoints, scale):

    """
    This function ...
    :param parameter_range:
    :param npoints:
    :param scale:
    :return: 
    """

    # Generate the grid points, based on the scale
    if scale == "linear": values = parameter_range.linear(npoints, as_list=True)
    elif scale == "logarithmic": values = parameter_range.log(npoints, as_list=True)
    elif scale == "sqrt": values = parameter_range.sqrt(npoints, as_list=True)
    else: raise ValueError("Invalid scale: " + str(scale))

    # Re-arrange the list to contain the grid points from the center towards the edges
    values = sequences.rearrange_from_middle(values)

    # Return the values
    return values

# -----------------------------------------------------------------
