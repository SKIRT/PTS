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

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(GridModelGenerator, self).setup(**kwargs)

        # Prompt for most sampled parameters
        if self.most_sampled_parameters is None: self.prompt_most_sampled_parameters()

        # Prompt for sampling weights
        if self.sampling_weights is None: self.prompt_sampling_weights()

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

    @property
    def parameter_scales(self):

        """
        This function ...
        :return: 
        """

        scales = dict()

        # Loop over the free parameter
        for label in self.fitting_run.free_parameter_labels:

            # Check whether scales were given as input
            if self.scales is not None and label in self.scales: scales[label] = self.scales[label]
            elif self.config.scales is not None and label in self.config.scales: scales[label] = self.config.scales[label]
            else: raise ValueError("Scale was not set for " + label)

        # Return the scales
        return scales

    # -----------------------------------------------------------------

    def scale_for_parameter(self, label):

        """
        This function ...
        :param label: 
        :return: 
        """

        # Check whether scales were given as input
        if self.scales is not None and label in self.scales: return self.scales[label]
        elif self.config.scales is not None and label in self.config.scales: return self.config.scales[label]
        else: raise ValueError("Scale was not set for " + label)

    # -----------------------------------------------------------------

    def generate(self):

        """
        THis function ...
        :return: 
        """

        # Inform the user
        log.info("Generating the model grid points ...")

        # Generate the grid points (as dictionary of lists)
        grid_points = self.generate_grid_points_different_scales(self.parameter_scales, most_sampled=self.config.most_sampled_parameters, weights=self.config.sampling_weights) # returns dictionary

        # Convert into lists
        grid_points = self.grid_points_to_lists(grid_points)

        # Create iterator of combinations
        iterator = sequences.iterate_lists_combinations(grid_points)

        # Create name iterator
        name_iterator = strings.alphabet_strings_iterator()

        # Generate the initial parameter sets
        # Loop over the number of required models minus the number of fixed model parameter sets
        for index in range(self.config.nmodels):

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
