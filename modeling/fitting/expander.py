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

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.basics.log import log
from ...core.tools.utils import lazyproperty, memoize_method
from ...core.tools import sequences
from ...core.tools import nr, numbers
from ...core.basics.containers import DefaultOrderedDict

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

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 6. Generate the model parameters
        self.generate_models()

        # 7. Set the paths to the input files
        #if self.needs_input: self.set_input()

        # 8. Adjust the ski template
        #self.adjust_ski()

        # 9. Fill the tables for the current generation
        #self.fill_tables()

        # Launch the models
        self.launch()

        # Show
        self.show()

        # Write
        self.write()

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

    @property
    def initial_parameter_values(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.first_guess_parameter_values

    # -----------------------------------------------------------------

    @lazyproperty
    def generation(self):

        """
        Thisfunction ...
        :return:
        """

        return self.fitting_run.get_generation(self.config.generation)

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
        for label in self.parameter_labels:
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

    @property
    def chi_squared_table(self):

        """
        This function ...
        :return:
        """

        return self.generation.chi_squared_table

    # -----------------------------------------------------------------

    @property
    def up(self):

        """
        This function ...
        :return:
        """

        return self.config.direction == up

    # -----------------------------------------------------------------

    @property
    def down(self):

        """
        Thisn function ...
        :return:
        """

        return self.config.direction == down

    # -----------------------------------------------------------------

    @property
    def both(self):

        """
        This function ...
        :return:
        """

        return self.config.direction == both

    # -----------------------------------------------------------------

    @lazyproperty
    def series(self):

        """
        This function ...
        :return:
        """

        if self.up: return numbers.get_linear_series(self.config.npoints, start=1, step=1)
        elif self.down: return numbers.get_linear_series(self.config.npoints, start=-1, step=-1)
        elif self.both: return numbers.get_alternating_series(self.config.npoints, start=1)
        else: raise ValueError("Invalid direction")

    # -----------------------------------------------------------------

    def generate_models(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the new model parameters ...")

        # Loop over the parameters
        for label in self.parameter_labels:

            # Get sorted unique values
            unique_values = sequences.ordered(self.unique_parameter_values_scalar[label])
            lowest_value = unique_values[0]
            highest_value = unique_values[-1]

            # Get step
            if self.is_linear(label):

                steps = numbers.differences(unique_values)
                step = sequences.get_all_close_value(steps)

                # Add new values
                #for i in range(self.config.npoints):
                #    if self.up:
                #    elif self.down:
                #    self.new_values.append(highest_value + i * step)

                if self.up:

                    for i in self.series:
                        new_value = highest_value + i * step
                        self.new_parameter_values[label].append(new_value)

                elif self.down:

                    for i in self.series:
                        new_value = lowest_value + i * step
                        self.new_parameter_values[label].append(new_value)

                elif self.both:

                    for i, value in zip(self.series, sequences.alternate([highest_value, lowest_value], self.config.npoints)):
                        new_value = value + i * step
                        self.new_parameter_values[label].append(new_value)

                else: raise ValueError("Invalid direction")

            elif self.is_logarithmic(label):

                factors = numbers.quotients(unique_values)
                factor = sequences.get_all_close_value(factors)
                print(factor)



            else: raise ValueError("Invalid scale")

    # -----------------------------------------------------------------

    def launch(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching ...")

        # Create manager
        #self.manager = SimulationManager()

        # Run the manager
        #self.manager.run(assignment=assignment, timing=self.timing_table, memory=self.memory_table)

        # status=status, info_tables=[parameters, chi_squared], remotes=remotes, simulations=simulations)

        # Set the actual number of simulations for this generation
        #self.generation_info.nsimulations = self.nmodels

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

# -----------------------------------------------------------------
