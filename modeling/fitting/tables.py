#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.generations Contains the GenerationsTable class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from ...core.basics.table import SmartTable
from ...core.tools import tables

# -----------------------------------------------------------------

class BestParametersTable(SmartTable):

    """
    This class ...
    """

    column_info = [("Generation name", str, None, "Name for the generation")]

    # -----------------------------------------------------------------

    @classmethod
    def initialize(cls, parameters):

        """
        This function ...
        :param parameters:
        :return:
        """

        # Set the labels
        for label in parameters: cls.column_info.append((label, float, None, "value for " + label))

        # Call the initialize function of the best parameters table function
        return super(BestParametersTable, cls).initialize()

    # -----------------------------------------------------------------

    @property
    def parameter_labels(self):

        """
        This function ...
        :return:
        """

        labels = []

        for name in self.colnames:
            if name == "Generation name": continue
            labels.append(name)

        return labels

    # -----------------------------------------------------------------

    def add_entry(self, generation_name, parameter_values):

        """
        This function ...
        :param generation_name:
        :param parameter_values:
        :return:
        """

        values = [generation_name]

        # Add the parameter values
        for label in self.parameter_labels: values.append(parameter_values[label])

        # Resize string columns for the new values
        self._resize_string_columns(values)

        # Add a row to the table
        self.add_row(values)

# -----------------------------------------------------------------

class GenerationsTable(SmartTable):

    """
    This class ...
    """

    column_info = [("Generation name", str, None, "Name for the generation"),
                   ("Generation index", int, None, "Index of the generation"),
                   ("Launching time", str, None, "Time of launching the generation simulations"),
                   ("Method", str, None, "Method used for model generation"),
                   ("Wavelength grid level", int, None, "level of the wavelength gid"),
                   ("Dust grid level", int, None, "level of the dust grid"),
                   ("Number of simulations", int, None, "number of simulations (individuals) in the generation"),
                   ("Self-absorption", bool, None, "dust self-absorption enabled")]

    # -----------------------------------------------------------------

    @classmethod
    def initialize(cls, parameters):

        """
        This function ...
        :param parameters:
        :return:
        """

        for label in parameters:

            cls.column_info.append(("Minimum value for " + label, float, None, "minimum value for " + label))
            cls.column_info.append(("Maximum value for " + label, float, None, "Maximum value for " + label))

        cls.column_info.append(("Finishing time", str, None, "Time of finishing the generation"))

        # Call the initialize function of the generations table function
        return super(GenerationsTable, cls).initialize()

    # -----------------------------------------------------------------

    @property
    def parameter_labels(self):

        """
        This function ...
        :return:
        """

        labels = []

        for name in self.colnames:
            if name.startswith("Minimum value for "):
                labels.append(name.split("value for ")[1])

        return labels

    # -----------------------------------------------------------------

    @property
    def generation_names(self):

        """
        This function ...
        :return:
        """

        return list(self["Generation name"])

    # -----------------------------------------------------------------

    def set_finishing_time(self, generation_name, timestamp):

        """
        This function ...
        :param generation_name:
        :param timestamp:
        :return:
        """

        if generation_name not in self.generation_names: raise ValueError("Generation '" + generation_name + "' does not exist in the table")

        for i in range(len(self)):

            if self["Generation name"][i] == generation_name: # match

                self._resize_string_column("Finishing time", timestamp)

                # Set the value
                self["Finishing time"].mask[i] = False
                self["Finishing time"][i] = timestamp

                break

    # -----------------------------------------------------------------

    def is_finished(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        # Find the index in the table for the specified generation
        index = tables.find_index(self, generation_name, "Generation name")

        # Return whether the value for the finishing time of the generation is not masked
        return not self["Finishing time"].mask[index]

    # -----------------------------------------------------------------

    def add_entry(self, name, index, timestamp, method, wavelength_grid_level, dust_grid_level, nsimulations, selfabsorption, ranges):

        """
        This function ...
        :param name:
        :param index:
        :param timestamp:
        :param method:
        :param wavelength_grid_level:
        :param dust_grid_level:
        :param nsimulations:
        :param selfabsorption:
        :param ranges:
        :return:
        """

        values = [name, index, timestamp, method, wavelength_grid_level, dust_grid_level, nsimulations, selfabsorption]

        # Add the boundaries (min, max) of the parameter ranges as seperate column values
        for label in self.parameter_labels:
            values.append(ranges[label].min)
            values.append(ranges[label].max)

        # Add None for the finishing time
        values.append(None)

        # Resize string columns for the new values
        self._resize_string_columns(values)

        # Add a row to the table
        self.add_row(values)

# -----------------------------------------------------------------

class ParametersTable(SmartTable):

    """
    This class ...
    """

    column_info = [("Simulation name", str, None, "name of the simulation")]

    # -----------------------------------------------------------------

    @classmethod
    def initialize(cls, parameters):

        """
        This function ...
        :param parameters:
        :return:
        """

        for label in parameters:
            cls.column_info.append((label, float, None, "value for " + label))

        # Call the initialize function of the parameters table function
        return super(ParametersTable, cls).initialize()

    # -----------------------------------------------------------------

    @property
    def parameter_labels(self):

        """
        This function ...
        :return:
        """

        labels = []

        for name in self.colnames:
            if name == "Simulation name": continue
            labels.append(name)

        return labels

    # -----------------------------------------------------------------

    def parameter_values_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        values = dict()

        index = tables.find_index(self, simulation_name, "Simulation name")

        for name in self.colnames:
            if name == "Simulation name": continue
            values[name] = self[name][index]

        # Return the values
        return values

    # -----------------------------------------------------------------

    def add_entry(self, name, parameter_values):

        """
        This function ...
        :param name:
        :param parameter_values:
        :return:
        """

        values = [name]

        # Add the parameter values
        for label in self.parameter_labels: values.append(parameter_values[label])

        # Resize string columns for the new values
        self._resize_string_columns(values)

        # Add a row to the table
        self.add_row(values)

# -----------------------------------------------------------------

class ChiSquaredTable(SmartTable):

    """
    This class ...
    """

    column_info = [("Simulation name", str, None, "name of the simulation"),
                   ("Chi squared", float, None, "chi-squared value")]

    # -----------------------------------------------------------------

    @property
    def best_simulation_name(self):

        """
        This function ...
        :return:
        """

        index = np.argmin(self["Chi squared"])
        return self["Simulation name"][index]

    # -----------------------------------------------------------------

    @property
    def best_simulation_name_and_chi_squared(self):

        """
        This function ...
        :return:
        """

        index = np.argmin(self["Chi squared"])
        return self["Simulation name"][index], self["Chi squared"][index]

    # -----------------------------------------------------------------

    def chi_squared_for(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        index = tables.find_index(self, simulation_name, "Simulation name")
        return self["Chi squared"][index]

    # -----------------------------------------------------------------

    def add_entry(self, name, chi_squared):

        """
        This function ...
        :param name:
        :param chi_squared:
        :return:
        """

        values = [name, chi_squared]

        # Resize string columns for the new values
        self._resize_string_columns(values)

        # Add a row to the table
        self.add_row(values)

# -----------------------------------------------------------------
