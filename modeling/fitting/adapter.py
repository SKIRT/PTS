#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.adapter Contains the FittingAdapter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.basics.log import log
from ...core.basics.configuration import prompt_mapping, prompt_variable, prompt_automatic
from .configuration import magnitude_ranges, percentual_ranges, default_magnitude_range, default_percentual_range, get_initial_value
from ...core.basics.range import QuantityRange
from ...core.tools import sequences

# -----------------------------------------------------------------

class FittingAdapter(FittingComponent):
    
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
        super(FittingAdapter, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The fitting run
        self.fitting_run = None

        # The configuration
        self.configuration = None

    # -----------------------------------------------------------------

    @property
    def representation(self):

        """
        This function ...
        :return:
        """

        return "representation" in self.config.properties

    # -----------------------------------------------------------------

    @property
    def filters(self):

        """
        This function ...
        :return:
        """

        return "filters" in self.config.properties

    # -----------------------------------------------------------------

    @property
    def ranges(self):

        """
        This function ...
        :return:
        """

        return "ranges" in self.config.properties

    # -----------------------------------------------------------------

    @property
    def genetic(self):

        """
        This function ...
        :return:
        """

        return "genetic" in self.config.properties

    # -----------------------------------------------------------------

    @property
    def grid(self):

        """
        This function ...
        :return:
        """

        return "grid" in self.config.properties

    # -----------------------------------------------------------------

    @property
    def units(self):

        """
        This function ...
        :return:
        """

        return "units" in self.config.properties

    # -----------------------------------------------------------------

    @property
    def types(self):

        """
        This function ...
        :return:
        """

        return "types" in self.config.properties

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Representation
        if self.representation: self.adapt_representation()

        # Filters
        if self.filters: self.adapt_filters()

        # Ranges
        if self.ranges: self.adapt_ranges()

        # Genetic
        if self.genetic and self.genetic_fitting: self.adapt_genetic()

        # Grid
        if self.grid and self.grid_fitting: self.adapt_grid()

        # Units
        if self.units: self.adapt_units()

        # Types
        if self.types: self.adapt_types()

        # 13. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(FittingAdapter, self).setup()

        # Load the fitting run
        self.fitting_run = self.load_fitting_run(self.config.name)

        # Load the configuration
        self.configuration = self.fitting_run.fitting_configuration

    # -----------------------------------------------------------------

    @property
    def free_parameter_labels(self):

        """
        This function ...
        :return:
        """

        return self.configuration.free_parameters

    # -----------------------------------------------------------------

    @property
    def fitting_method(self):

        """
        This function ...
        :return:
        """

        return self.configuration.method

    # -----------------------------------------------------------------

    @property
    def grid_fitting(self):

        """
        This function ...
        :return:
        """

        return self.fitting_method == "grid"

    # -----------------------------------------------------------------

    @property
    def genetic_fitting(self):

        """
        Thisn function ...
        :return:
        """

        return self.fitting_method == "genetic"

    # -----------------------------------------------------------------

    @property
    def model_name(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.model_name

    # -----------------------------------------------------------------

    @property
    def model_definition(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.model_definition

    # -----------------------------------------------------------------

    @property
    def model_representation_names(self):

        """
        This function ...
        :return:
        """

        return self.static_model_suite.representations_for_model(self.model_name)

    # -----------------------------------------------------------------

    def adapt_representation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adapting initial model representation ...")

        # Prompt
        name = prompt_automatic("initial_representation", "initial model representation", default=self.configuration.initial_representation, choices=self.model_representation_names)

        # Changed?
        changed = name != self.configuration.initial_representation
        if changed and self.config.save:
            self.configuration.initial_representation = name
            self.configuration.save()

    # -----------------------------------------------------------------

    def adapt_filters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adapting fitting filters ...")

        # Prompt
        filter_names = prompt_variable("filters", "string_list", "fitting filters", default=self.configuration.filters)

        # Changed?
        changed = not sequences.same_contents(filter_names, self.configuration.filters)
        if changed and self.config.save:
            self.configuration.filters = filter_names
            self.configuration.save()

    # -----------------------------------------------------------------

    def adapt_ranges(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adapting parameter ranges ...")

        # Loop over the free parameters
        for label in self.free_parameter_labels:

            # Get the initial value for this parameter
            value = get_initial_value(self.model_definition, label)

            # Use magnitudes
            if label in magnitude_ranges:

                # Get default magnitude
                magnitude = magnitude_ranges[label]

                # Construct the actual range
                suggestion = QuantityRange.within_magnitude(value, magnitude)

            # Use percentage
            elif label in percentual_ranges:

                # Get default
                fraction = float(percentual_ranges[label]) / 100.

                # Construct the actual range
                rel_min = 1. - fraction
                rel_max = 1. + fraction
                suggestion = QuantityRange.around(value, rel_min=rel_min, rel_max=rel_max)

            # Use default percentage
            else:

                # Get default
                fraction = float(default_percentual_range) / 100.

                # Construct the actual range
                rel_min = 1. - fraction
                rel_max = 1. + fraction
                suggestion = QuantityRange.around(value, rel_min=rel_min, rel_max=rel_max)

            # Set suggestions
            suggestions = [suggestion]

            # Get new range
            property_name = label + "_range"
            default_range = self.configuration[property_name]
            parameter_range = prompt_automatic(property_name, "range for '" + label +  "'", default=default_range, suggestions=suggestions)

            # Changed?
            changed = parameter_range != default_range
            if changed and self.config.save:
                self.configuration[property_name] = parameter_range
                self.configuration.save()

    # -----------------------------------------------------------------

    def adapt_genetic(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adapting genetic settings ...")

        # Adapt
        changed = prompt_mapping(self.configuration["genetic"], contains=self.config.contains, not_contains=self.config.not_contains, exact_name=self.config.exact_name, exact_not_name=self.config.exact_not_name, startswith=self.config.startswith, endswith=self.config.endswith)

        # Save
        if changed and self.config.save: self.configuration.save()

    # -----------------------------------------------------------------

    def adapt_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adapting grid settings ...")

        # Adapt
        changed = prompt_mapping(self.configuration["grid"], contains=self.config.contains, not_contains=self.config.not_contains, exact_name=self.config.exact_name, exact_not_name=self.config.exact_not_name, startswith=self.config.startswith, endswith=self.config.endswith)

        # Save
        if changed and self.config.save: self.configuration.save()

    # -----------------------------------------------------------------

    def adapt_units(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adapting units ...")

        # Adapt
        changed = prompt_mapping(self.configuration["units"])

        # Save
        if changed and self.config.save: self.configuration.save()

    # -----------------------------------------------------------------

    def adapt_types(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adapting types ...")

        # Adapt
        changed = prompt_mapping(self.configuration["types"])

        # Save
        if changed and self.config.save: self.configuration.save()

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

# -----------------------------------------------------------------
