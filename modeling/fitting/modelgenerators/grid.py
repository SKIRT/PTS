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
from ....core.basics.log import log
from .generator import ModelGenerator
from ....core.tools import sequences
from ....core.tools import strings
from ....core.basics.configuration import prompt_string_list, prompt_weights

# -----------------------------------------------------------------

class GridModelGenerator(ModelGenerator):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(GridModelGenerator, self).__init__(*args, **kwargs)

        # Number of sampling points per parameter
        self.npoints = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(GridModelGenerator, self).setup(**kwargs)

        # Get npoints
        if "npoints" in kwargs: self.npoints = kwargs.pop("npoints")

        # Prompt for most sampled parameters
        if not self.has_npoints and self.most_sampled_parameters is None: self.prompt_most_sampled_parameters()

        # Prompt for sampling weights
        if not self.has_npoints and self.sampling_weights is None: self.prompt_sampling_weights()

    # -----------------------------------------------------------------

    @property
    def has_npoints(self):

        """
        This function ...
        :return:
        """

        return self.npoints is not None

    # -----------------------------------------------------------------

    def prompt_most_sampled_parameters(self):

        """
        THis function ...
        :return: 
        """

        # Inform the user
        log.info("Prompting for the most sampled parameters ...")

        # Get the most sampled parameters
        self.most_sampled_parameters = prompt_string_list("parameters", "free parameter(s) which get the most sampling points", choices=self.fitting_run.free_parameter_labels, required=False)

    # -----------------------------------------------------------------

    def prompt_sampling_weights(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Prompting for the sampling weights ...")

        # Get the sampling weights
        self.sampling_weights = prompt_weights("weights", "relative sampling for the free parameters: " + ",".join(self.fitting_run.free_parameter_labels), required=False)

    # -----------------------------------------------------------------

    # @property
    # def parameter_scales(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     scales = dict()
    #
    #     # Loop over the free parameter
    #     for label in self.fitting_run.free_parameter_labels:
    #
    #         # Check whether scales were given as input
    #         if self.scales is not None and label in self.scales: scales[label] = self.scales[label]
    #         elif self.config.scales is not None and label in self.config.scales: scales[label] = self.config.scales[label]
    #         else: #raise ValueError("Scale was not set for '" + label + "'")
    #             # Take from grid fitting configuration
    #             scales[label] = self.fitting_run.grid_settings[label + "_scale"]
    #
    #     # Return the scales
    #     return scales

    # -----------------------------------------------------------------

    def scale_for_parameter(self, label):

        """
        This function ...
        :param label: 
        :return: 
        """

        # Check whether scales were given as input
        #if self.scales is not None and label in self.scales: return self.scales[label]
        #elif self.config.scales is not None and label in self.config.scales: return self.config.scales[label]
        #else: raise ValueError("Scale was not set for '" + label + "'")
        return self.scales[label]

    # -----------------------------------------------------------------

    @property
    def combined_npoints(self):

        """
        This function ...
        :return:
        """

        if not self.has_npoints: raise ValueError("Number of points per parameter not defined")
        total = 1
        for label in self.npoints:
            npoints = self.npoints[label]
            total *= npoints
        return total

    # -----------------------------------------------------------------

    @property
    def target_nmodels(self):

        """
        This function ...
        :return:
        """

        if self.has_npoints: return self.combined_npoints
        else: return self.config.nmodels

    # -----------------------------------------------------------------

    def generate(self):

        """
        THis function ...
        :return: 
        """

        # Inform the user
        log.info("Generating the model parameters ...")

        # Generate the grid points (as dictionary of lists)
        grid_points = self.generate_grid_points()
        #print(grid_points)

        # Create iterator of combinations
        iterator = sequences.iterate_lists_combinations(*grid_points)

        # Create name iterator
        name_iterator = strings.alphabet_strings_iterator()

        # Generate the initial parameter sets
        # Loop over the number of required models minus the number of fixed model parameter sets
        #print("nmodels:", self.target_nmodels)
        for index in range(self.target_nmodels):

            # The next combination
            parameters_model = list(iterator.next())  # returns tuple

            # Generate a new individual name
            name = name_iterator.next()

            # Loop over all the parameters
            for i in range(len(parameters_model)):

                # Get the parameter value
                value = parameters_model[i]

                # Get the parameter label
                label = self.fitting_run.free_parameter_labels[i]

                # Add the parameter value to the dictionary
                self.parameters[label][name] = value

    # -----------------------------------------------------------------

    def generate_grid_points(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Generating grid points ...")

        # Generate as dictionary
        points_per_parameter = self.generate_grid_points_different_scales(self.scales, npoints=self.npoints,
                                                                         most_sampled=self.config.most_sampled_parameters,
                                                                         weights=self.config.sampling_weights)  # returns dictionary

        # Convert into lists
        return self.grid_points_to_lists(points_per_parameter)

    # -----------------------------------------------------------------

    @property
    def individual_names(self):

        """
        This function ...
        :return: 
        """

        one_free_parameter = self.fitting_run.free_parameter_labels[0]
        return self.parameters[one_free_parameter].keys()

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Writing ...")

# -----------------------------------------------------------------
