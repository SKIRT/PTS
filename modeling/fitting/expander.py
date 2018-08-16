#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.expander Contains the ParameterExpander class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.basics.log import log
from ...core.tools.utils import lazyproperty, memoize_method
from ...core.tools import sequences
from ...core.tools import nr, numbers, strings, types
from ...core.basics.containers import DefaultOrderedDict
from ...core.tools import formatting as fmt
from ...core.tools.stringify import tostr, stringify_dict
from .evaluate import get_parameter_values_for_named_individual
from .manager import GenerationManager
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

up = "up"
down = "down"
both = "both"
directions = [up, down, both]

# -----------------------------------------------------------------

class ParameterExpander(FittingComponent):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(ParameterExpander, self).__init__(*args, **kwargs)

        # The generation info
        self.info = None

        # The new parameter values
        self.new_parameter_values = DefaultOrderedDict(list)

        # The model parameters
        self.parameters = DefaultOrderedDict(OrderedDict)

        # The new simulation names
        self.new_simulation_names = []

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Generate the parameter values
        self.generate_parameters()

        # Generate the models
        self.generate_models()

        # Fill the tables for the current generation
        self.fill_tables()

        # Show
        self.show()

        # Write
        self.write()

        # Launch the models
        self.launch()

        # Update the tables
        self.update()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ParameterExpander, self).setup(**kwargs)

        # Load the generation info
        self.info = self.fitting_run.get_generation_info(self.config.generation)

    # -----------------------------------------------------------------

    @property
    def nnew_simulations(self):

        """
        This function ...
        :return:
        """

        return len(self.new_simulation_names)

    # -----------------------------------------------------------------

    @property
    def ncurrent_simulations(self):

        """
        This function ...
        :return:
        """

        return self.generation.nsimulation_directories

    # -----------------------------------------------------------------

    @property
    def ntotal_simulations(self):

        """
        This function ...
        :return:
        """

        return self.ncurrent_simulations + self.nnew_simulations

    # -----------------------------------------------------------------

    @property
    def fitting_run_name(self):

        """
        This function ...
        :return:
        """

        return self.config.run

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_run(self):

        """
        This function ...
        :return:
        """

        return self.load_fitting_run(self.config.run)

    # -----------------------------------------------------------------

    @property
    def free_parameter_labels(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.free_parameter_labels

    # -----------------------------------------------------------------

    @property
    def nfree_parameters(self):

        """
        This function ...
        :return:
        """

        return len(self.free_parameter_labels)

    # -----------------------------------------------------------------

    @property
    def parameter_ranges(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.free_parameter_ranges

    # -----------------------------------------------------------------

    @lazyproperty
    def parameter_labels(self):

        """
        This function ...
        :return:
        """

        if self.config.parameters is not None: return self.config.parameters
        else: return self.free_parameter_labels

    # -----------------------------------------------------------------

    @property
    def nparameters(self):

        """
        This function ...
        :return:
        """

        return len(self.parameter_labels)

    # -----------------------------------------------------------------

    @property
    def has_multiple_parameters(self):

        """
        This function ...
        :return:
        """

        return self.nparameters > 1

    # -----------------------------------------------------------------

    @property
    def has_one_parameter(self):

        """
        This function ...
        :return:
        """

        return self.nparameters == 1

    # -----------------------------------------------------------------

    @property
    def has_all_parameters(self):

        """
        This function ...
        :return:
        """

        return self.nparameters == self.nfree_parameters

    # -----------------------------------------------------------------

    @property
    def parameter_units(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.parameter_units

    # -----------------------------------------------------------------

    def get_parameter_unit(self, label):

        """
        Thisf unction ...
        :param label:
        :return:
        """

        return self.parameter_units[label]

    # -----------------------------------------------------------------

    def has_parameter_unit(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        return label in self.parameter_units and self.parameter_units[label] is not None

    # -----------------------------------------------------------------

    @property
    def initial_parameter_values(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.first_guess_parameter_values

    # -----------------------------------------------------------------

    @property
    def generation_name(self):

        """
        This function ...
        :return:
        """

        return self.config.generation

    # -----------------------------------------------------------------

    @lazyproperty
    def generation(self):

        """
        Thisfunction ...
        :return:
        """

        return self.fitting_run.get_generation(self.generation_name)

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
    def grid_settings(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.grid_settings

    # -----------------------------------------------------------------

    @lazyproperty
    def parameter_scales(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary for the scales
        scales = dict()

        # Get the scales for each free parameter
        for label in self.free_parameter_labels:
            key = label + "_scale"
            scales[label] = self.grid_settings[key]

        # Return the scales dict
        return scales

    # -----------------------------------------------------------------

    def is_linear(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        return self.parameter_scales[label] == "linear"

    # -----------------------------------------------------------------

    def is_logarithmic(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        return self.parameter_scales[label] == "logarithmic"

    # -----------------------------------------------------------------

    @property
    def parameters_table(self):

        """
        This function ...
        :return:
        """

        return self.generation.parameters_table

    # -----------------------------------------------------------------

    @lazyproperty
    def simulation_names(self):

        """
        This function ...
        :return:
        """

        return self.parameters_table.simulation_names

    # -----------------------------------------------------------------

    @property
    def first_simulation_name(self):

        """
        This function ...
        :return:
        """

        return self.simulation_names[0]

    # -----------------------------------------------------------------

    @lazyproperty
    def simulation_indices(self):

        """
        This function ...
        :return:
        """

        indices = []
        for name in self.simulation_names:
            if "__" not in name: continue
            index = int(name.split("__")[-1])
            indices.append(index)
        return indices

    # -----------------------------------------------------------------

    @property
    def has_simulation_indices(self):

        """
        This function ...
        :return:
        """

        return len(self.simulation_indices) > 0

    # -----------------------------------------------------------------

    @lazyproperty
    def last_simulation_index(self):

        """
        This function ...
        :return:
        """

        if self.has_simulation_indices: return max(self.simulation_indices)
        else: return -1

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_parameter_values(self):

        """
        This function ...
        :return:
        """

        return self.parameters_table.unique_parameter_values

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_parameter_values_scalar(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary
        values_scalar = DefaultOrderedDict(list)

        # Loop over the parameters
        for label in self.unique_parameter_values:
            for value in self.unique_parameter_values[label]:
                scalar_value = value.to(self.get_parameter_unit(label)).value
                values_scalar[label].append(scalar_value)

        # Return the scalar values
        return values_scalar

    # -----------------------------------------------------------------

    @memoize_method
    def get_nunique_parameter_values(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        return len(self.unique_parameter_values[label])

    # -----------------------------------------------------------------

    @memoize_method
    def get_sorted_unique_parameter_values_scalar(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        # Get sorted unique values
        return sequences.ordered(self.unique_parameter_values_scalar[label])

    # -----------------------------------------------------------------

    @memoize_method
    def get_lowest_unique_parameter_value_scalar(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        return self.get_sorted_unique_parameter_values_scalar(label)[0]

    # -----------------------------------------------------------------

    @memoize_method
    def get_lowest_unique_parameter_value(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        return self.get_lowest_unique_parameter_value_scalar(label) * self.get_parameter_unit(label)

    # -----------------------------------------------------------------

    @memoize_method
    def get_highest_unique_parameter_value_scalar(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        return self.get_sorted_unique_parameter_values_scalar(label)[-1]

    # -----------------------------------------------------------------

    @memoize_method
    def get_highest_unique_parameter_value(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        return self.get_highest_unique_parameter_value_scalar(label) * self.get_parameter_unit(label)

    # -----------------------------------------------------------------

    @property
    def individuals_table(self):

        """
        This function ...
        :return:
        """

        return self.generation.individuals_table

    # -----------------------------------------------------------------

    @lazyproperty
    def individual_names(self):

        """
        This function ...
        :return:
        """

        return self.individuals_table.individual_names

    # -----------------------------------------------------------------

    @property
    def chi_squared_table(self):

        """
        This function ...
        :return:
        """

        return self.generation.chi_squared_table

    # -----------------------------------------------------------------

    @memoize_method
    def up(self, label):

        """
        This function ...
        :return:
        """

        if types.is_string_type(self.config.direction): return self.config.direction == up
        elif types.is_dictionary(self.config.direction): return self.config.direction[label] == up
        else: raise ValueError("Invalid type for 'direction'")

    # -----------------------------------------------------------------

    @memoize_method
    def down(self, label):

        """
        Thisn function ...
        :return:
        """

        if types.is_string_type(self.config.direction): return self.config.direction == down
        elif types.is_dictionary(self.config.direction): return self.config.direction[label] == down
        else: raise ValueError("Invalid type for 'direction'")

    # -----------------------------------------------------------------

    @memoize_method
    def both(self, label):

        """
        This function ...
        :return:
        """

        if types.is_string_type(self.config.direction): return self.config.direction == both
        elif types.is_dictionary(self.config.direction): return self.config.direction[label] == both
        else: raise ValueError("Invalid type for 'direction'")

    # -----------------------------------------------------------------

    @memoize_method
    def npoints(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        if types.is_integer_type(self.config.npoints): return self.config.npoints
        elif types.is_dictionary(self.config.npoints): return self.config.npoints[label]
        else: raise ValueError("Invalid type for 'npoints'")

    # -----------------------------------------------------------------

    @memoize_method
    def series(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        if self.up(label): return numbers.get_linear_series(self.npoints(label), start=1, step=1)
        elif self.down(label): return numbers.get_linear_series(self.npoints(label), start=-1, step=-1)
        elif self.both(label): return numbers.get_alternating_series(self.npoints(label), start=1)
        else: raise ValueError("Invalid direction")

    # -----------------------------------------------------------------

    def add_new_parameter_value(self, label, value):

        """
        This function ...
        :param label:
        :param value:
        :return:
        """

        # Debugging
        log.debug("Adding new parameter value of " + tostr(value) + " for parameter '" + label + "' ...")

        # Add unit
        value = value * self.get_parameter_unit(label)

        # Add
        self.new_parameter_values[label].append(value)

    # -----------------------------------------------------------------

    def generate_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the new parameters values ...")

        # Loop over the parameters
        for label in self.parameter_labels:

            # Linear range
            if self.is_linear(label): self.generate_parameters_linear(label)

            # Logarithmic range
            elif self.is_logarithmic(label): self.generate_parameters_logarithmic(label)

            # Invalid scale
            else: raise ValueError("Invalid scale")

    # -----------------------------------------------------------------

    @memoize_method
    def get_grid_step(self, label):

        """
        This function ...
        :param label: 
        :return: 
        """

        # Get the unique values
        unique_values = self.get_sorted_unique_parameter_values_scalar(label)

        # Determine the step size
        steps = numbers.differences(unique_values)
        return sequences.get_all_close_value(steps)

    # -----------------------------------------------------------------

    def generate_parameters_linear(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        # Inform the user
        log.info("Generating new parameter values for '" + label + "' on a linear scale ...")

        # Determine the step size
        step = self.get_grid_step(label)

        # Show
        log.debug("The step size of the linear grid of '" + label + "' is " + str(step) + " " + tostr(self.get_parameter_unit(label)))

        # Get lowest and highest value
        lowest_value = self.get_lowest_unique_parameter_value_scalar(label)
        highest_value = self.get_highest_unique_parameter_value_scalar(label)

        # Debug
        log.debug("Current lowest value: " + tostr(lowest_value))
        log.debug("Current highest value: " + tostr(highest_value))

        # Add points above the current range
        if self.up(label):

            for i in self.series(label):
                new_value = highest_value + i * step
                self.add_new_parameter_value(label, new_value)

        # Add points below the current range
        elif self.down(label):

            for i in self.series(label):
                new_value = lowest_value + i * step
                self.add_new_parameter_value(label, new_value)

        # Add points both above and below the current range
        elif self.both(label):

            for i, value in zip(self.series(label), sequences.alternate([highest_value, lowest_value], self.config.npoints)):
                new_value = value + i * step
                self.add_new_parameter_value(label, new_value)

        # Invalid direction
        else: raise ValueError("Invalid direction")

    # -----------------------------------------------------------------

    @memoize_method
    def get_grid_factor(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        # Get unique values
        unique_values = self.get_sorted_unique_parameter_values_scalar(label)

        # Determine the factor
        factors = numbers.quotients(unique_values)
        return sequences.get_all_close_value(factors)

    # -----------------------------------------------------------------

    def generate_parameters_logarithmic(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        # Inform the user
        log.info("Generating new parameter values for '" + label + "' on a logarithmic scale ...")

        # Get the grid factor
        factor = self.get_grid_factor(label)

        # Show
        log.debug("The increment of the logarithmic grid of '" + label + "' is a factor of " + str(factor))

        # Get lowest and highest value
        lowest_value = self.get_lowest_unique_parameter_value_scalar(label)
        highest_value = self.get_highest_unique_parameter_value_scalar(label)

        # Debug
        log.debug("Current lowest value: " + tostr(lowest_value))
        log.debug("Current highest value: " + tostr(highest_value))

        # Add points above the current range
        if self.up(label):

            for i in self.series(label):
                new_value = highest_value * factor ** i
                self.add_new_parameter_value(label, new_value)

        # Add points below the current range
        elif self.down(label):

            for i in self.series(label):
                new_value = lowest_value * factor ** i
                self.add_new_parameter_value(label, new_value)

        # Add points both above and below the current range
        elif self.both(label):

            for i, value in zip(self.series(label), sequences.alternate([highest_value, lowest_value], self.config.npoints)):
                new_value = value * factor ** i
                self.add_new_parameter_value(label, new_value)

        # Invalid
        else: raise ValueError("Invalid direction")

    # -----------------------------------------------------------------

    @lazyproperty
    def name_iterator(self):

        """
        This function ...
        :return:
        """

        # Create name iterator, increment
        name_iterator = strings.alphabet_strings_iterator()
        name_iterator.increment_to(self.individual_names)

        # Return
        return name_iterator

    # -----------------------------------------------------------------

    def generate_models(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the new models ...")

        # Generate models with one new parameter value combined with original parameter values
        self.generate_models_single_new_parameter()

        # Generate models with multiple new parameter values combined with original parameter values
        if self.has_multiple_parameters and not self.nfree_parameters == 2: self.generate_models_multiple_new_parameters()

        # Generate models with all new parameter values
        if self.has_multiple_parameters and self.has_all_parameters: self.generate_models_all_new_parameters()

    # -----------------------------------------------------------------

    def generate_models_single_new_parameter(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating models with a single new parameter and other original parameter values ...")

        # Loop over the parameters with new values
        for label in self.parameter_labels: self.generate_new_models_for_parameter(label)

    # -----------------------------------------------------------------

    @memoize_method
    def get_other_parameters(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        # Get labels
        if types.is_string_type(label): labels = [label]
        elif types.is_string_sequence(label) or types.is_string_tuple(label): labels = label
        else: raise ValueError("Invalid input")

        # Return other
        return sequences.get_other(self.free_parameter_labels, labels)

    # -----------------------------------------------------------------

    @memoize_method
    def get_grid_points_dict_for_parameters(self, parameter_labels):

        """
        This function ...
        :param parameter_labels:
        :return:
        """

        # Initialize
        grid_points = OrderedDict()

        # Loop over the parameters
        for parameter_label in parameter_labels: grid_points[parameter_label] = self.get_sorted_unique_parameter_values_scalar(parameter_label)

        # Return
        return grid_points

    # -----------------------------------------------------------------

    @memoize_method
    def get_new_parameter_values(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        return self.new_parameter_values[label]

    # -----------------------------------------------------------------

    @memoize_method
    def get_new_parameter_values_scalar(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        return [value.to(self.get_parameter_unit(label)).value for value in self.new_parameter_values[label]]

    # -----------------------------------------------------------------

    @memoize_method
    def get_new_grid_points_dict_for_parameters(self, parameter_labels):

        """
        This function ...
        :param parameter_labels:
        :return:
        """

        # Initialize
        grid_points = OrderedDict()

        # Loop over the parameters
        for parameter_label in parameter_labels: grid_points[parameter_label] = self.get_new_parameter_values_scalar(parameter_label)

        # Return
        return grid_points

    # -----------------------------------------------------------------

    @memoize_method
    def get_grid_points_dict_for_other_parameters(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        # Initialize
        grid_points = OrderedDict()

        # Loop over the other parameters
        for parameter_label in self.get_other_parameters(label): grid_points[parameter_label] = self.get_sorted_unique_parameter_values_scalar(parameter_label)

        # Return
        return grid_points

    # -----------------------------------------------------------------

    def generate_new_models_for_parameter(self, label):

        """
        This function ...
        :param label: 
        :return: 
        """

        # Debugging
        log.debug("Generating models for the new parameter values of '" + label + "' ...")

        # Loop over the new parameter values
        for value in self.new_parameter_values[label]:

            # Debugging
            log.debug("Generating new models with " + label + " = " + tostr(value) + " ...")

            # Other parameter values
            parameters = self.get_grid_points_dict_for_other_parameters(label)
            other_labels = parameters.keys()
            #print(parameter_labels)

            # Get number of models
            nmodels = 1
            for other_label in other_labels: nmodels *= self.get_nunique_parameter_values(other_label)

            # Debugging
            log.debug("Number of models to be generated: " + str(nmodels))

            # Create iterator of combinations
            iterator = sequences.iterate_lists_combinations(*parameters.values())

            # Loop over the grid points of the other parameters
            for index in range(nmodels):

                # The next combination
                other_values = list(iterator.next())  # returns tuple

                # Generate a new individual name
                name = self.name_iterator.next()

                # Add the parameter value to the dictionary
                self.parameters[label][name] = value.to(self.get_parameter_unit(label)).value
                for other_label, other_value in zip(other_labels, other_values): self.parameters[other_label][name] = other_value

    # -----------------------------------------------------------------

    def get_nnew_parameter_values(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        return len(self.new_parameter_values[label])

    # -----------------------------------------------------------------

    def generate_models_multiple_new_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating models with multiple new parameter values combined with old parameter values ...")

        # Set lengths (combine at least 2 parameters, and at most all of the parameters, unless these are ALL free parameters; then the function below will already take care of that)
        lengths = list(range(2,min(self.nparameters, self.nfree_parameters-1)+1))

        # Debugging
        log.debug("Combining the parameters '" + tostr(self.parameter_labels) + "' in sets of " + tostr(lengths) + " ...")

        # Loop over the combinations of multiple
        for parameter_labels in sequences.combinations(self.parameter_labels, lengths):

            # Convert into tuple
            parameter_labels = tuple(parameter_labels)

            # Debugging
            log.debug("Generating models with the new values for the parameters " + tostr(parameter_labels) + " ...")

            # Get new parameter values
            new_parameters = self.get_new_grid_points_dict_for_parameters(parameter_labels)

            # Other parameter values
            other_parameters = self.get_grid_points_dict_for_other_parameters(parameter_labels)
            other_labels = other_parameters.keys()

            # Get number of models
            nnew_combinations = 1
            for label in parameter_labels: nnew_combinations *= self.get_nnew_parameter_values(label)
            nother_combinations = 1
            for other_label in other_labels: nother_combinations *= self.get_nunique_parameter_values(other_label)
            nmodels = nnew_combinations * nother_combinations

            # Debugging
            log.debug("Number of models to be generated for this combination: " + str(nmodels) + " (" + str(nnew_combinations) + " combinations of new parameters and " + str(nother_combinations) + " combinations of other parameters)")

            # Create iterator of combinations of new parameter values
            new_iterator = sequences.iterate_lists_combinations(*new_parameters.values())

            # Loop over the grid points of the new parameters
            for i in range(nnew_combinations):

                # The next combination
                new_values = list(new_iterator.next())  # returns tuple

                # Create iterator of combinations of other parameter values
                other_iterator = sequences.iterate_lists_combinations(*other_parameters.values())

                # Loop over the grid points of the other parameters
                for j in range(nother_combinations):

                    # The next combination
                    other_values = list(other_iterator.next())  # returns tuple

                    # Generate a new individual name
                    name = self.name_iterator.next()

                    # Add the parameter value to the dictionary
                    for new_label, new_value in zip(parameter_labels, new_values): self.parameters[new_label][name] = new_value
                    for other_label, other_value in zip(other_labels, other_values): self.parameters[other_label][name] = other_value

    # -----------------------------------------------------------------

    def generate_models_all_new_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating models with all new parameter values ...")

        # Get number of models
        nmodels = 1
        for label in self.free_parameter_labels: nmodels *= self.get_nnew_parameter_values(label)

        # Debugging
        log.debug("Number of models to be generated: " + str(nmodels))

        # Get new parameter values
        new_parameters = self.get_new_grid_points_dict_for_parameters(tuple(self.free_parameter_labels))
        new_parameter_labels = new_parameters.keys()

        # Create iterator of combinations of new parameter values
        iterator = sequences.iterate_lists_combinations(*new_parameters.values())

        # Loop over the new models
        for i in range(nmodels):

            # The next combination
            values = list(iterator.next())  # returns tuple

            # Generate a new individual name
            name = self.name_iterator.next()

            # Loop over the parameters
            for parameter_label, value in zip(new_parameter_labels, values):

                # Add the parameter value to the dictionary
                self.parameters[parameter_label][name] = value

    # -----------------------------------------------------------------

    @lazyproperty
    def new_individual_names(self):

        """
        This function ...
        :return:
        """

        return self.parameters[self.parameters.keys()[0]].keys()

    # -----------------------------------------------------------------

    @property
    def model_names(self):

        """
        This function ...
        :return:
        """

        return self.new_individual_names

    # -----------------------------------------------------------------

    def fill_tables(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Filling the tables for the current generation ...")

        # Loop over the model names
        counter = 0
        for name in self.model_names:

            # Get the simulation index
            simulation_index = self.last_simulation_index + 1 + counter

            # Generate the simulation name
            simulation_name = self.object_name + "__" + self.fitting_run_name + "__" + self.generation_name + "__" + str(simulation_index)

            # Add the new simulation name
            self.new_simulation_names.append(simulation_name)

            # Debugging
            log.debug("Adding an entry to the individuals table with:")
            log.debug("")
            log.debug(" - Simulation name: " + simulation_name)
            log.debug(" - Individual_name: " + name)
            log.debug("")

            # Add entry
            self.individuals_table.add_entry(simulation_name, name)

            # Get the parameter values
            parameter_values = get_parameter_values_for_named_individual(self.parameters, name, self.fitting_run)

            # Debugging
            log.debug("Adding entry to the parameters table with:")
            log.debug("")
            log.debug(" - Simulation name: " + simulation_name)
            for label in parameter_values: log.debug(" - " + label + ": " + tostr(parameter_values[label], scientific=True, fancy=True, ndigits=self.fitting_run.ndigits_dict[label]))
            log.debug("")

            # Add an entry to the parameters table
            self.parameters_table.add_entry(simulation_name, parameter_values)

            # Increment counter
            counter += 1

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

        # Show
        self.show_parameter_values()

        # Show the models
        self.show_models()

    # -----------------------------------------------------------------

    def has_new_parameter_values(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        return label in self.new_parameter_values and len(self.new_parameter_values[label]) > 0

    # -----------------------------------------------------------------

    @memoize_method
    def get_new_parameter_values_above(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        # No new values for this parameter
        if not self.has_new_parameter_values(label): return []

        values = []
        highest = self.get_highest_unique_parameter_value(label)
        for value in self.new_parameter_values[label]:
            if value > highest: values.append(value)
        return values

    # -----------------------------------------------------------------

    @memoize_method
    def get_new_parameter_values_below(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        # No new values for this parameter
        if not self.has_new_parameter_values(label): return []

        values = []
        lowest = self.get_lowest_unique_parameter_value(label)
        for value in self.new_parameter_values[label]:
            if value < lowest: values.append(value)
        return values

    # -----------------------------------------------------------------

    def show_parameter_values(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the parameter values ...")

        # Loop over the free parameters
        for label in self.free_parameter_labels:

            # Get the unit
            unit = self.get_parameter_unit(label)

            print("")
            print(fmt.underlined + fmt.green + label + fmt.reset + " [" + tostr(unit) + "]:")
            print(self.parameter_scales[label])
            print("")

            # Show values below original range
            below = self.get_new_parameter_values_below(label)
            for value in below: print("  " + fmt.cyan + tostr(value.to(unit).value) + fmt.reset)

            # Show original values
            for value in self.get_sorted_unique_parameter_values_scalar(label): print("  " + tostr(value))

            # Show values above original range
            above = self.get_new_parameter_values_above(label)
            for value in above: print("  " + fmt.cyan + tostr(value.to(unit).value) + fmt.reset)

        print("")

    # -----------------------------------------------------------------

    def show_models(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the model parameters ...")

        # Print in columns
        with fmt.print_in_columns() as print_row:

            # Set column names and units
            column_names = ["Individual"] + self.free_parameter_labels
            column_units = [""] + [self.get_parameter_unit(label) for label in self.free_parameter_labels]

            # Show the header
            print_row(*column_names)
            if not sequences.all_none(column_units): print_row(*column_units)

            # Loop over the new individuals
            for name in self.new_individual_names:

                row = []
                row.append(name)

                # Add parameter value
                for label in self.free_parameter_labels: row.append(self.parameters[label][name])

                # Show row
                print_row(*row)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # 2. Write the individuals table
        self.write_individuals()

        # 3. Write the parameters table
        self.write_parameters()

    # -----------------------------------------------------------------

    def write_individuals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the individuals table ...")

        # Determine the path
        test_path = fs.join(self.generation_path, "individuals_new.dat")

        # Save
        self.individuals_table.saveto(test_path)

    # -----------------------------------------------------------------

    def write_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the parameters table ...")

        # Determine the path
        test_path = fs.join(self.generation_path, "parameters_new.dat")

        # Save
        self.parameters_table.saveto(test_path)

    # -----------------------------------------------------------------

    def launch(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching the new simulations ...")

        # Create simulation manager
        manager = GenerationManager()

        # Set fitting run and generation name
        manager.config.run = self.fitting_run_name
        manager.config.generation = self.generation_name

        # Set options
        #manager.config.extra = config.extra
        #manager.config.offline = config.offline
        #manager.config.lazy = config.lazy
        #manager.config.find_simulations = config.find_simulations
        #manager.config.find_remotes = config.find_remotes
        #manager.config.produce_missing = config.produce_missing
        #manager.config.check_paths = config.check_paths
        #manager.config.correct_paths = config.correct_paths
        #manager.config.confirm_correction = config.confirm_correction
        #manager.config.fix_success = config.fix_success
        #manager.config.check_analysis = config.check_analysis
        #manager.config.write_status = config.write_status
        #manager.config.correct_status = config.correct_status

        # Not interactive
        manager.config.interactive = False

        # Initialize empty list of commands
        manager.config.commands = []

        # EXTRA
        #extra = " --parallelization 16:2:2 --not_mimic_scheduling --scheduling/walltime 54000 --scheduling/nodes 1"
        #extra = " --parallelization 64:8:1" # For nancy

        # Set extra options
        extra = " --not_mimic_scheduling"
        if self.config.host is not None: extra += " --host '" + self.config.host.as_string() + "'"
        if self.config.parallelization is not None: extra += " --parallelization '" + self.config.parallelization.as_string() + "'"
        if self.config.walltime is not None: extra += " --scheduling/walltime " + str(self.config.walltime) # in seconds
        if self.config.nnodes is not None: extra += " --scheduling/nodes " + str(self.config.nnodes)

        # Set the mimic commands
        for simulation_name in self.new_simulation_names:

            # Get the parameter values
            parameter_values = self.parameters_table.parameter_values_for_simulation(simulation_name)
            parameter_values_string = stringify_dict(parameter_values, quote_character="'", identity_symbol=":")[1]

            # Construct mimic command
            mimic_command = 'mimic "' + self.first_simulation_name + '" "' + simulation_name + '" ' + '--labeled "' + parameter_values_string + '"' + extra
            #print(mimic_command)

            # Debugging
            log.debug("Command for running simulation '" + simulation_name + "': '" + mimic_command + "'")

            # Add command
            manager.config.commands.append(mimic_command)

        # Run the generation manager
        manager.run()
        if manager.has_failed: raise RuntimeError("Something went wrong launching the simulations")

    # -----------------------------------------------------------------

    def update(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Updating the generation tables ...")

        # Individuals table
        if self.config.update_individuals: self.update_individuals()

        # Parameters table
        if self.config.update_parameters: self.update_parameters()

        # Generation info
        if self.config.update_info: self.update_info()

        # Generations table
        if self.config.update_generations: self.update_generations()

    # -----------------------------------------------------------------

    def update_individuals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Updating the individuals table ...")

        # Determine the path
        filepath = self.generation.individuals_table_path
        new_filepath = fs.join(self.generation_path, "individuals_new.dat")

        # Backup the old individuals table
        fs.backup_file(filepath, suffix="old", exists="backup")
        fs.remove_file(filepath)

        # Save
        self.individuals_table.saveto(filepath)

        # Remove new, it has now become unnecessary
        if fs.is_file(new_filepath): fs.remove_file(new_filepath)

    # -----------------------------------------------------------------

    def update_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Updating the parameters table ...")

        # Determine the path
        filepath = self.generation.parameters_table_path
        new_filepath = fs.join(self.generation_path, "parameters_new.dat")

        # Backup the old parameters table
        fs.backup_file(filepath, suffix="old", exists="backup")
        fs.remove_file(filepath)

        # Save
        self.parameters_table.saveto(filepath)

        # Remove the new, it has now become unnecessary
        if fs.is_file(new_filepath): fs.remove_file(new_filepath)

    # -----------------------------------------------------------------

    def update_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Updating the generation info ...")

        # Set the new number of models
        self.info.nsimulations = self.ntotal_simulations

        # Save the info
        self.info.save()

    # -----------------------------------------------------------------

    @property
    def generations_table(self):
        return self.fitting_run.generations_table

    # -----------------------------------------------------------------

    def update_generations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Updating the generations table ...")

        # Set the new number of simulations
        self.generations_table.set_nsimulations(self.generation_name, self.ntotal_simulations)

        # TODO: update the ranges!!

        # Save the table
        self.generations_table.save()

# -----------------------------------------------------------------
