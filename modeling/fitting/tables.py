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
from ...core.basics.range import RealRange

# -----------------------------------------------------------------

class BestParametersTable(SmartTable):

    """
    This class ...
    """

    column_info = [("Generation name", str, None, "Name for the generation")]

    # -----------------------------------------------------------------

    @classmethod
    def initialize(cls, parameters, units):

        """
        This function ...
        :param parameters:
        :param units:
        :return:
        """

        # Set the labels
        for label in parameters:
            unit = units[label] if label in units else None
            cls.column_info.append((label, float, unit, "value for " + label))

        # Add the Chi squared column
        cls.column_info.append(("Chi squared", float, None, "chi-squared value"))

        # Call the initialize function of the best parameters table function
        return super(BestParametersTable, cls).initialize()

    # -----------------------------------------------------------------

    @property
    def generation_names(self):

        """
        This function ...
        :return:
        """

        return list(self["Generation name"])

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

    def add_entry(self, generation_name, parameter_values, chi_squared):

        """
        This function ...
        :param generation_name:
        :param parameter_values:
        :param chi_squared:
        :return:
        """

        values = [generation_name]

        # Add the parameter values
        for label in self.parameter_labels: values.append(parameter_values[label])

        # Add the chi squared value
        values.append(chi_squared)

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
                   ("Number of photon packages", int, None, "number of photon packages per wavelength"),
                   ("Self-absorption", bool, None, "dust self-absorption enabled"),
                   ("Transient heating", bool, None, "transient heating enabled")]

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

    @property
    def genetic_generations(self):

        """
        This function ...
        :return:
        """

        # Initialize a list to contain the generation names
        names = []

        # Loop over the rows in the table
        for i in range(len(self)):

            # Add the generation name
            if not self["Generation index"].mask[i]: names.append(self["Generation name"][i])

        # Return the generation names
        return names

    # -----------------------------------------------------------------

    @property
    def genetic_generations_with_initial(self):

        """
        This function ...
        :return:
        """

        # Get the list of genetic generation names without the initial generation
        names = self.genetic_generations

        # Prepend the name of the initial generation if present
        if "initial" in self["Generation name"]: names = ["initial"] + names

        # Return the list of generation name
        return names

    # -----------------------------------------------------------------

    @property
    def finished_generations(self):

        """
        This function ...
        :return:
        """

        # Initialize list
        names = []

        # Loop over each row (generation)
        for i in range(len(self)):

            name = self["Generation name"][i]
            finished = not self["Finishing time"].mask[i]

            if finished: names.append(name)

        # Return the list of finished generations
        return names

    # -----------------------------------------------------------------

    def parameter_ranges_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        # Find the row index of the specified generation
        index = tables.find_index(self, generation_name, "Generation name")

        # Create a dictionary to contain the ranges
        ranges = dict()

        # Loop over the parameter labels
        for label in self.parameter_labels:

            # Get the minimum and maximum value
            min_value = self["Minimum value for " + label][index]
            max_value = self["Maximum value for " + label][index]

            # Create the range
            range = RealRange(min_value, max_value)

            # Add the range to the dictionary
            ranges[label] = range

        # Return the ranges dictionary
        return ranges

    # -----------------------------------------------------------------

    def set_finishing_time(self, generation_name, timestamp):

        """
        This function ...
        :param generation_name:
        :param timestamp:
        :return:
        """

        # Check if the generation exists
        if generation_name not in self.generation_names: raise ValueError("Generation '" + generation_name + "' does not exist in the table")

        # Loop over the rows, find the entry for the specified generation
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

    def add_entry(self, name, index, timestamp, method, wavelength_grid_level, dust_grid_level, nsimulations, npackages, selfabsorption, transientheating, ranges):

        """
        This function ...
        :param name:
        :param index:
        :param timestamp:
        :param method:
        :param wavelength_grid_level:
        :param dust_grid_level:
        :param nsimulations:
        :param npackages:
        :param selfabsorption:
        :param transientheating:
        :param ranges:
        :return:
        """

        values = [name, index, timestamp, method, wavelength_grid_level, dust_grid_level, nsimulations, npackages, selfabsorption, transientheating]

        # Add the boundaries (min, max) of the parameter ranges as seperate column values
        for label in self.parameter_labels:
            values.append(ranges[label].min)
            values.append(ranges[label].max)

        # Add None for the finishing time
        values.append(None)

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
    def initialize(cls, parameters, units):

        """
        This function ...
        :param parameters:
        :param units:
        :return:
        """

        for label in parameters:
            unit = units[label] if label in units else None
            cls.column_info.append((label, float, unit, "value for " + label))

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

    @property
    def simulation_names(self):

        """
        This function ...
        :return:
        """

        return list(self["Simulation name"])

    # -----------------------------------------------------------------

    def add_entry(self, name, chi_squared):

        """
        This function ...
        :param name:
        :param chi_squared:
        :return:
        """

        # Set the values
        values = [name, chi_squared]

        # Add a row to the table
        self.add_row(values)

# -----------------------------------------------------------------

class ModelProbabilitiesTable(SmartTable):

    """
    This class ...
    """

    column_info = [("Simulation name", str, None, "name of the simulation")]

    # -----------------------------------------------------------------

    @property
    def parameter_labels(self):

        """
        This function ...
        :return:
        """

        labels = []

        for name in self.colnames:

            if name == "Simulation name" or name == "Probability": continue
            labels.append(name)

        return labels

    # -----------------------------------------------------------------

    @classmethod
    def initialize(cls, parameters, units):

        """
        This function ...
        :param parameters:
        :param units:
        :return:
        """

        # Add columns for the parameter values
        for label in parameters:

            unit = units[label] if label in units else None
            cls.column_info.append((label, float, unit, "value for " + label))

        # Add column for probabilities
        cls.column_info.append(("Probability", float, None, "model probability"))

        # Call the initialize function of the parameters table function
        return super(ModelProbabilitiesTable, cls).initialize()

    # -----------------------------------------------------------------

    def add_entry(self, simulation_name, parameter_values, probability):

        """
        This function ...
        :param simulation_name:
        :param probability:
        :return:
        """

        # Set the values
        values = [simulation_name]

        # Add the parameter values
        for label in self.parameter_labels:
            values.append(parameter_values[label])

        # Add the probability
        values.append(probability)

        # Add a row to the table
        self.add_row(values)

# -----------------------------------------------------------------

class ParameterProbabilitiesTable(SmartTable):

    """
    This class ...
    """

    column_info = [("Parameter value", float, None, "value of the parameter"),
                   ("Probability", float, None, "probability for this parameter value")]

    # -----------------------------------------------------------------

    def add_entry(self, value, probability):

        """
        This function ...
        :param value:
        :param probability:
        :return:
        """

        # Set the values
        values = [value, probability]

        # Add a row to the table
        self.add_row(values)

# -----------------------------------------------------------------
