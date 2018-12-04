#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.tables Contains table classes: WeightsTable, BestParametersTable, GenerationsTable,
#  ParametersTable, ChiSquaredTable, ModelProbabilitiesTable, and ParameterProbabilitiesTable

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ...core.basics.table import SmartTable
from ...core.tools import tables
from ...core.basics.range import RealRange, QuantityRange
from ...core.basics.curve import FilterCurve
from ...core.units.parsing import parse_unit as u
from ...core.tools import time
from ...core.tools import sequences
from ...core.basics.distribution import Distribution
from ...core.tools import numbers
from ...core.tools.utils import lazyproperty

# -----------------------------------------------------------------

class RunsTable(SmartTable):

    """
    This class ...
    """

    # Add column info
    _column_info = OrderedDict()
    _column_info["Run name"] = (str, None, "Name for the fitting run")
    _column_info["Model name"] = (str, None, "Name of the model used")

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(RunsTable, self).__init__(*args, **kwargs)

        # Add column info
        self.add_all_column_info(self._column_info)

    # -----------------------------------------------------------------

    def add_run(self, run):

        """
        This function ...
        :param run:
        :return:
        """

        self.add_entry(run.name, run.model_name)

    # -----------------------------------------------------------------

    def add_entry(self, run_name, model_name):

        """
        This function ...
        :param run_name:
        :param model_name:
        :return:
        """

        values = [run_name, model_name]

        # Add row
        self.add_row(values)

    # -----------------------------------------------------------------

    @property
    def run_names(self):
        return list(self["Run name"])

    # -----------------------------------------------------------------

    def model_for_run(self, run_name):

        """
        This function ...
        :param run_name:
        :return:
        """

        index = tables.find_index(self, run_name)
        if index is None: raise ValueError("Fitting run '" + run_name + "' not found in the table")
        return self["Model name"][index]

# -----------------------------------------------------------------

class WeightsTable(FilterCurve):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Set kwargs
        kwargs["y_name"] = "Weight"
        kwargs["y_description"] = "Weight given to the filter"
        kwargs["y_unit"] = None

        # Call the constructor of the base class
        super(WeightsTable, self).__init__(*args, **kwargs)

# -----------------------------------------------------------------

class BestParametersTable(SmartTable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Check
        if "parameters" in kwargs: from_astropy = False
        else: from_astropy = True

        # Get properties
        if not from_astropy:
            parameters = kwargs.pop("parameters")
            units = kwargs.pop("units")
        else: parameters = units = None

        # Call the constructor of the base class
        super(BestParametersTable, self).__init__(*args, **kwargs)

        # Not from constructor called by Astropy
        if not from_astropy:

            # Generation column
            self.add_column_info("Generation name", str, None, "Name for the generation")

            # Parameters columns
            for label in parameters:
                unit = units[label] if label in units else None
                self.add_column_info(label, float, unit, "value for " + label)

            # Add the chi squared column
            self.add_column_info("Chi squared", float, None, "chi-squared value")

    # -----------------------------------------------------------------

    @property
    def generation_names(self):
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
            if name == "Generation name" or name == "Chi squared": continue
            labels.append(name)

        return labels

    # -----------------------------------------------------------------

    def has_generation(self, generation_name):

        """
        This function ...
        :param generation_name: 
        :return: 
        """

        return generation_name in self.generation_names

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

    def remove_entry(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        index = tables.find_index(self, generation_name, "Generation name")
        self.remove_row(index)

# -----------------------------------------------------------------

class GenerationsTable(SmartTable):

    """
    This class ...
    """

    _column_info = OrderedDict()
    _column_info["Generation name"] = (str, None, "name for the generation")
    _column_info["Generation index"] = (int, None, "index of the generation")
    _column_info["Launching time"] = (str, None, "time of launching the generation simulations")
    _column_info["Method"] = (str, None, "method used for model generation")
    _column_info["Wavelength grid name"] = (str, None, "name of the wavelength gid")
    _column_info["Model representation"] = (str, None, "representation of the model")
    _column_info["Number of simulations"] = (int, None, "number of simulations (individuals) in the generation")
    _column_info["Number of photon packages"] = (int, None, "number of photon packages per wavelength")
    _column_info["Self-absorption"] = (bool, None, "dust self-absorption enabled")
    _column_info["Transient heating"] = (bool, None, "transient heating enabled")
    _column_info["Spectral convolution"] = (bool, None, "spectral convolution enabled")
    _column_info["Use images"] = (bool, None, "use images")

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Check
        if "parameters" in kwargs: from_astropy = False
        else: from_astropy = True

        # Get properties
        if not from_astropy:
            parameters = kwargs.pop("parameters")
            units = kwargs.pop("units")
        else: parameters = units = None

        # Call the constructor of the base class
        super(GenerationsTable, self).__init__(*args, **kwargs)

        # If not called from astropy
        if not from_astropy:

            # Add column info
            self.add_all_column_info(self._column_info)

            # RANGES
            # Loop over the parameters
            for label in parameters:

                # Get the unit
                unit = u(units[label]) if label in units and units[label] is not None else None

                # Add columns for the minimum and maximum value for this parameter
                self.add_column_info("Minimum value for " + label, float, unit, "Minimum value for " + label)
                self.add_column_info("Maximum value for " + label, float, unit, "Maximum value for " + label)

            # SCALES
            # Loop over the parameters
            for label in parameters: self.add_column_info("Scale for " + label, str, None, "Scale for " + label) # Add column for the scale of this parameter

            # Add finishing time column
            self.add_column_info("Finishing time", str, None, "Time of finishing the generation")

    # -----------------------------------------------------------------

    def index_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        row_index = tables.find_index(self, generation_name)
        return self["Generation index"][row_index]

    # -----------------------------------------------------------------

    @property
    def ngenerations(self):

        """
        This function ...
        :return:
        """

        return len(self)

    # -----------------------------------------------------------------

    @property
    def parameter_labels(self):

        """
        This function ...
        :return:
        """

        if len(self.colnames) == 0: self._setup()

        labels = []

        for name in self.colnames:
            if name.startswith("Minimum value for "):
                labels.append(name.split("value for ")[1])

        return labels

    # -----------------------------------------------------------------

    @property
    def generation_names(self):
        return list(self["Generation name"])

    # -----------------------------------------------------------------

    @property
    def last_generation_name(self):
        if len(self) > 0: return self["Generation name"][-1]
        else: return None

    # -----------------------------------------------------------------

    @property
    def grid_generations(self):

        """
        This function ...
        :return:
        """

        # Initialize a list to contain the generation names
        names = []

        # Loop over the rows in the table
        for i in range(len(self)):

            # Add the generation name
            if "grid" in self["Generation name"][i]: names.append(self["Generation name"][i])

        # Return the generation names
        return names

    # -----------------------------------------------------------------

    @property
    def ngrid_generations(self):
        return len(self.grid_generations)

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
    def ngenetic_generations(self):
        return len(self.genetic_generations)

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
    def genetic_generations_with_initial_names_and_indices(self):

        """
        This function ...
        :return:
        """

        names = self.genetic_generations_with_initial

        names_and_indices = []

        for name in names:

            if name == "initial": index = -1
            else: index = self.index_for_generation(name)

            names_and_indices.append((name, index))

        return names_and_indices

    # -----------------------------------------------------------------

    @property
    def ngenetic_generations_with_initial(self):
        return len(self.genetic_generations_with_initial)

    # -----------------------------------------------------------------

    @property
    def all_finished(self):

        """
        This function ...
        :return:
        """

        # Loop over each row (generation)
        for i in range(len(self)):

            finished = not self["Finishing time"].mask[i]

            # At least one is not finished
            if not finished: return False

        # All finished
        return True

    # -----------------------------------------------------------------

    @property
    def has_finished(self):

        """
        This function ...
        :return:
        """

        # Loop over each row (generation)
        for i in range(len(self)):

            finished = not self["Finishing time"].mask[i]

            # At least one is finished
            if finished: return True

        # None are finished
        return False

    # -----------------------------------------------------------------

    @property
    def has_unfinished(self):

        """
        This function ...
        :return:
        """

        # Loop over each generation
        for i in range(len(self)):

            finished = not self["Finishing time"].mask[i]

            # If one unfinished is encountered, return True
            if not finished: return True

        # None are unfinished
        return False

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

        # Create a dictionary to contain the ranges
        ranges = dict()

        # Loop over the parameter labels
        for label in self.parameter_labels:

            # Get the range
            range = self.parameter_range_for_generation(generation_name, label)

            # Add the range to the dictionary
            ranges[label] = range

        # Return the ranges dictionary
        return ranges

    # -----------------------------------------------------------------

    def parameter_range_for_generation(self, generation_name, label):

        """
        This function ...
        :param generation_name: 
        :param label: 
        :return: 
        """

        # Find the row index of the specified generation
        index = tables.find_index(self, generation_name, "Generation name")

        # Get the minimum and maximum value
        min_value = self["Minimum value for " + label][index]
        max_value = self["Maximum value for " + label][index]

        # Create the range
        column_unit = self.column_unit("Minimum value for " + label)
        if column_unit is None: range = RealRange(min_value, max_value)
        else: range = QuantityRange(min_value, max_value, unit=column_unit)

        # Return the range
        return range

    # -----------------------------------------------------------------

    def _get_generation_index(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        # COMPLETELY DIFFERENT FROM INDEX_FOR_GENERATION!!
        return self.find_index(generation_name, "Generation name")

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
        #for i in range(len(self)):
        index = self._get_generation_index(generation_name)

        # Resize column if necessary
        self._resize_string_column("Finishing time", timestamp)

        # Set the value
        self["Finishing time"].mask[index] = False
        self["Finishing time"][index] = timestamp

    # -----------------------------------------------------------------

    def set_nsimulations(self, generation_name, nsimulations):

        """
        This function ...
        :param generation_name:
        :param nsimulations:
        :return:
        """

        # Check if the generation exists
        if generation_name not in self.generation_names: raise ValueError("Generation '" + generation_name + "' does not exist in the table")

        # Get the index
        index = self._get_generation_index(generation_name)

        # Set the number of simulations
        self["Number of simulations"][index] = nsimulations

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

    def remove_last_entry(self):

        """
        This function ...
        :return:
        """

        self.remove_row(len(self)-1)

    # -----------------------------------------------------------------

    def remove_entry(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        index = tables.find_index(self, generation_name, "Generation name")
        self.remove_row(index)

    # -----------------------------------------------------------------

    def remove_entries(self, names):

        """
        This function ...
        :param names:
        :return:
        """

        indices = tables.find_indices(self, names, "Generation name")
        self.remove_rows(indices)

    # -----------------------------------------------------------------

    def remove_other_entries(self, names):

        """
        This function ...
        :param names:
        :return:
        """

        indices = tables.find_indices(self, names, "Generation name")
        other_indices = [index for index in range(len(self)) if index not in indices]
        self.remove_rows(other_indices)

    # -----------------------------------------------------------------

    def add_entry(self, generation_info, ranges, scales):

        """
        This function ...
        :param generation_info:
        :param ranges:
        :param scales:
        :return:
        """

        # Generate a timestamp
        timestamp = time.timestamp()

        # Add an entry to the generations table
        name = generation_info.name
        index = generation_info.index
        method = generation_info.method
        wg_name = generation_info.wavelength_grid_name
        representation = generation_info.model_representation_name
        nsimulations = generation_info.nsimulations
        npackages = generation_info.npackages
        selfabsorption = generation_info.selfabsorption
        transientheating = generation_info.transient_heating
        spectralconvolution = generation_info.spectral_convolution
        use_images = generation_info.use_images

        # Call other function
        self.add_entry_impl(name, index, timestamp, method, wg_name, representation, nsimulations, npackages, selfabsorption, transientheating, use_images, spectralconvolution, ranges, scales)

    # -----------------------------------------------------------------

    def add_entry_impl(self, name, index, timestamp, method, wavelength_grid_name, representation, nsimulations, npackages,
                  selfabsorption, transientheating, spectralconvolution, use_images, ranges, scales):

        """
        This function ...
        :param name:
        :param index:
        :param timestamp:
        :param method:
        :param wavelength_grid_name:
        :param representation:
        :param nsimulations:
        :param npackages:
        :param selfabsorption:
        :param transientheating:
        :param spectralconvolution:
        :param use_images:
        :param ranges:
        :param scales:
        :return:
        """

        values = [name, index, timestamp, method, wavelength_grid_name, representation, nsimulations, npackages, selfabsorption, transientheating, spectralconvolution, use_images]

        # Add the boundaries (min, max) of the parameter ranges as seperate column values
        for label in self.parameter_labels:
            values.append(ranges[label].min)
            values.append(ranges[label].max)

        # Add the scale
        for label in self.parameter_labels:
            values.append(scales[label])

        # Add None for the finishing time
        values.append(None)

        # Add a row to the table
        self.add_row(values)

# -----------------------------------------------------------------

class IndividualsTable(SmartTable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        THe constructor ...
        :param args: 
        :param kwargs: 
        """

        # Call the constructor of the base class
        super(IndividualsTable, self).__init__(*args, **kwargs)

        # Add column info
        self.add_column_info("Simulation name", str, None, "name of the simulation")
        self.add_column_info("Individual name", str, None, "name of the individual")

    # -----------------------------------------------------------------

    @property
    def simulation_names(self):
        return list(self["Simulation name"])

    # -----------------------------------------------------------------

    @property
    def individual_names(self):
        return list(self["Individual name"])

    # -----------------------------------------------------------------

    def has_individual(self, individual_name):

        """
        This function ...
        :param individual_name: 
        :return: 
        """

        index = tables.find_index(self, individual_name, "Individual name")
        return index is not None

    # -----------------------------------------------------------------

    def get_simulation_name(self, individual_name):

        """
        This function ...
        :param individual_name: 
        :return: 
        """

        index = tables.find_index(self, individual_name, "Individual name")
        if index is None: raise ValueError("Individual '" + individual_name + "' not found in this table")
        return self["Simulation name"][index]

    # -----------------------------------------------------------------

    def get_individual_name(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        index = tables.find_index(self, simulation_name)
        if index is None: raise ValueError("Simulation '" + simulation_name + "' not found in this table")
        return self["Individual name"][index]

    # -----------------------------------------------------------------

    def add_entry(self, simulation_name, individual_name):

        """
        This function ...
        :param simulation_name:
        :param individual_name:
        :return:
        """

        if len(self.colnames) == 0: self._setup()

        # Set the values
        values = [simulation_name, individual_name]

        # Add a row to the table
        self.add_row(values)

# -----------------------------------------------------------------

class ParametersTable(SmartTable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        if "parameters" in kwargs: from_astropy = False
        else: from_astropy = True

        if not from_astropy:
            parameters = kwargs.pop("parameters")
            units = kwargs.pop("units")
        else: parameters = units = None

        # Call the constructor of the base class
        super(ParametersTable, self).__init__(*args, **kwargs)

        # Set column info
        if not from_astropy:

            # Add simuation name column
            self.add_column_info("Simulation name", str, None, "name of the simulation")

            # Loop over the parameters
            for label in parameters:

                unit = units[label] if label in units else None
                self.add_column_info(label, float, unit, "value for " + label)

    # -----------------------------------------------------------------

    def unit_for(self, parameter_label):

        """
        This function ...
        :param parameter_label:
        :return:
        """

        return self.column_unit(parameter_label)

    # -----------------------------------------------------------------

    @property
    def simulation_names(self):
        return list(self["Simulation name"])

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

        return sequences.ordered(labels)

    # -----------------------------------------------------------------

    @property
    def parameter_units(self):
        units = dict()
        for label in self.parameter_labels: units[label] = self.unit_for(label)
        return units

    # -----------------------------------------------------------------

    def index_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Find index of row corresponding with the simulation
        index = tables.find_index(self, simulation_name, "Simulation name")
        if index is None: raise ValueError("Simulation not found in the table")
        return index

    # -----------------------------------------------------------------

    def get_simulation_name(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_value("Simulation name", index)

    # -----------------------------------------------------------------

    def parameter_values_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Get index
        index = self.index_for_simulation(simulation_name)

        # Return the parameter values
        return self.get_parameter_values(index)

    # -----------------------------------------------------------------

    def get_parameter_values(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        # Return the parameter values
        return self.get_values(self.parameter_labels, index, as_dict=True)

    # -----------------------------------------------------------------

    def simulation_for_parameter_values(self, values, return_index=False):

        """
        This function ...
        :param values:
        :param return_index:
        :return:
        """

        # Get list of values and list of parameter labels
        labels = values.keys()
        scalar_values = [values[label].to(self.column_unit(label)).value for label in values]

        # Find index of the simulation
        index = tables.find_one_index(self, scalar_values, labels)

        # Return the simulation name
        if return_index: return index
        else: return self.get_simulation_name(index)

    # -----------------------------------------------------------------

    def closest_simulation_for_parameter_values(self, values, return_index=False):

        """
        This function ...
        :param values:
        :param return_index:
        :return:
        """

        # Get list of values and list of parameter labels
        labels = values.keys()
        scalar_values = [values[label].to(self.column_unit(label)).value for label in values]

        # Find closest index of the simulation
        index = tables.find_closest_index(self, scalar_values, labels, rel_or_abs="rel")

        # Return the simulation name
        if return_index: return index
        else: return self.get_simulation_name(index)

    # -----------------------------------------------------------------

    def closest_simulations_for_parameter_values(self, values, nsimulations=2, return_indices=False):

        """
        This function ...
        :param values:
        :param nsimulations:
        :param return_indices:
        :return:
        """

        # Get list of values and list of parameter labels
        labels = values.keys()
        scalar_values = [values[label].to(self.column_unit(label)).value for label in values]

        # Find closest indices
        indices = tables.find_closest_indices(self, scalar_values, labels, nindices=nsimulations, rel_or_abs="rel")

        # Return
        if return_indices: return indices
        else: return [self.get_simulation_name(index) for index in indices]

    # -----------------------------------------------------------------

    @property
    def unique_parameter_values(self):

        """
        This function ...
        :return:
        """

        values = OrderedDict()

        for name in self.colnames:
            if name == "Simulation name": continue
            unique_values = sequences.unique_values(self[name])
            if self.column_unit(name) is not None: unique_values = [value * self.column_unit(name) for value in unique_values]
            values[name] = unique_values

        # Return
        return values

    # -----------------------------------------------------------------

    def add_entry(self, name, parameter_values):

        """
        This function ...
        :param name:
        :param parameter_values:
        :return:
        """

        if len(self.colnames) == 0: self._setup()

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

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ChiSquaredTable, self).__init__(*args, **kwargs)

        # Add column info
        self.add_column_info("Simulation name", str, None, "name of the simulation")
        self.add_column_info("Chi squared", float, None, "chi-squared value")

    # -----------------------------------------------------------------

    @property
    def worst_simulation_index(self):
        return np.argmax(self["Chi squared"])

    # -----------------------------------------------------------------

    @property
    def worst_simulation_name(self):
        return self["Simulation name"][self.worst_simulation_index]

    # -----------------------------------------------------------------

    def get_worst_simulation_names(self, nsimulations):

        """
        This function ...
        :param nsimulations:
        :return:
        """

        indices = np.argsort(self["Chi squared"])
        indices_reversed = indices[::-1]
        names = []
        for j in range(nsimulations):
            index = indices_reversed[j]
            name = self["Simulation name"][index]
            names.append(name)
        return names

    # -----------------------------------------------------------------

    @property
    def worst_chi_squared(self):
        return np.max(self["Chi squared"])

    # -----------------------------------------------------------------

    @property
    def best_simulation_index(self):
        return np.argmin(self["Chi squared"])

    # -----------------------------------------------------------------

    @property
    def best_simulation_name(self):
        return self["Simulation name"][self.best_simulation_index]

    # -----------------------------------------------------------------

    def get_best_simulation_names(self, nsimulations):

        """
        This function ...
        :param nsimulations:
        :return:
        """

        indices = np.argsort(self["Chi squared"])
        names = []
        for j in range(nsimulations):
            index = indices[j]
            name = self["Simulation name"][index]
            names.append(name)
        return names

    # -----------------------------------------------------------------

    @property
    def best_chi_squared(self):
        return np.min(self["Chi squared"])

    # -----------------------------------------------------------------

    @property
    def best_simulation_name_and_chi_squared(self):
        index = self.best_simulation_index
        return self["Simulation name"][index], self["Chi squared"][index]

    # -----------------------------------------------------------------

    def index_for_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return tables.find_index(self, simulation_name, "Simulation name")

    # -----------------------------------------------------------------

    def chi_squared_for(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        index = self.index_for_simulation(simulation_name)
        if index is None: return None
        return self["Chi squared"][index]

    # -----------------------------------------------------------------

    @property
    def simulation_names(self):
        return list(self["Simulation name"])

    # -----------------------------------------------------------------

    @property
    def chi_squared_values(self):
        return list(self["Chi squared"])

    # -----------------------------------------------------------------

    def has_simulation(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return name in self.simulation_names

    # -----------------------------------------------------------------

    def remove_simulation(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        index = self.index_for_simulation(name)
        self.remove_row(index)

    # -----------------------------------------------------------------

    def remove_all_simulations(self):

        """
        This function ...
        :return:
        """

        self.remove_all_rows()

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

    def set_chi_squared(self, simulation_name, chi_squared, index=None):

        """
        Thisn function ...
        :param simulation_name:
        :param chi_squared:
        :param index:
        :return:
        """

        # Get simulation index
        if index is None: index = self.find_index(simulation_name, "Simulation name")
        elif self.get_value("Simulation name", index) != simulation_name: raise ValueError("Invalid index for simulation '" + simulation_name + "'")

        # Set
        self.set_value("Chi squared", index, chi_squared)

    # -----------------------------------------------------------------

    def add_or_set_chi_squared(self, simulation_name, chi_squared):

        """
        This function ...
        :param simulation_name:
        :param chi_squared:
        :return:
        """

        index = self.find_index(simulation_name, "Simulation name")
        if index is None: self.add_entry(simulation_name, chi_squared)
        else: self.set_chi_squared(simulation_name, chi_squared, index=index)

    # -----------------------------------------------------------------

    def sort_as(self, simulation_names):

        """
        This function ...
        :param simulation_names: 
        :return: 
        """

        if len(simulation_names) != len(self): raise ValueError("Number of simulations does not have the same length as the table")
        if sequences.contains_duplicates(simulation_names): raise ValueError("Invalid list of simulation names: multiple occurences of the same name")

        new_column = []

        # Loop over the simulation names
        for name in simulation_names: new_column.append(self.chi_squared_for(name))

        # Replace the column
        self["Chi squared"] = new_column

    # -----------------------------------------------------------------

    @property
    def distribution(self):
        return Distribution.from_values("Chi squared", self["Chi squared"])

# -----------------------------------------------------------------

class ModelProbabilitiesTable(SmartTable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Check if from astropy
        if "parameters" in kwargs: from_astropy = False
        else: from_astropy = True

        # Get properties
        if not from_astropy:
            parameters = kwargs.pop("parameters")
            units = kwargs.pop("units")
        else: parameters = units = None

        # Call the constructor of the base class
        super(ModelProbabilitiesTable, self).__init__(*args, **kwargs)

        # Not from astropy
        if not from_astropy:

            # Add simulation name column
            self.add_column_info("Simulation name", str, None, "name of the simulation")

            # Loop over the parameters
            for label in parameters:

                unit = units[label] if label in units else None
                self.add_column_info(label, float, unit, "value for " + label)

            # Add column for the probabilities
            self.add_column_info("Probability", float, None, "model probability")

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

    def has_simulation(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        return simulation_name in self["Simulation name"]

    # -----------------------------------------------------------------

    @property
    def simulation_names(self):
        return list(self["Simulation name"])

    # -----------------------------------------------------------------

    @property
    def probabilities(self):
        return list(self["Probability"])

    # -----------------------------------------------------------------

    @property
    def nsimulations(self):
        return len(self["Simulation name"])

    # -----------------------------------------------------------------

    def add_entry(self, simulation_name, parameter_values, probability):

        """
        This function ...
        :param simulation_name:
        :param parameter_values:
        :param probability:
        :return:
        """

        if len(self.colnames) == 0: self._setup()

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

    @property
    def most_probable_simulation(self):

        """
        This function ...
        :return:
        """

        index = np.argmax(self["Probability"])
        return self.get_value("Simulation name", index)

    # -----------------------------------------------------------------

    @property
    def all_zero(self):
        return sequences.all_zero(self["Probability"])

    # -----------------------------------------------------------------

    @property
    def all_zero_but_one(self):
        return self.nzeros == self.nsimulations - 1

    # -----------------------------------------------------------------

    @property
    def has_zeros(self):
        return sequences.has_zero(self.probabilities)

    # -----------------------------------------------------------------

    @property
    def nzeros(self):
        return sequences.nzeros(self.probabilities)

# -----------------------------------------------------------------

class ParameterProbabilitiesTable(SmartTable):

    """
    This class ...
    """

    # Set column info
    _column_info = OrderedDict()
    _column_info["Value"] = (float, None, "value of the parameter")
    _column_info["Probability"] = (float, None, "probability for this parameter value")

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ParameterProbabilitiesTable, self).__init__(*args, **kwargs)

        # Add column info
        self.add_all_column_info(self._column_info)

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

    @property
    def most_probable_index(self):
        return np.argmax(self["Probability"])

    # -----------------------------------------------------------------

    @property
    def most_probable_value(self):
        return self["Value"][self.most_probable_index]

    # -----------------------------------------------------------------

    @property
    def values(self):
        return self.get_array("Value")

    # -----------------------------------------------------------------

    @property
    def probabilities(self):
        return self.get_array("Probability")

    # -----------------------------------------------------------------

    @lazyproperty
    def continuous_most_probable_value(self):
        return numbers.weighed_arithmetic_mean_numpy(self.values, self.probabilities)

    # -----------------------------------------------------------------

    @lazyproperty
    def most_probable_value_stddev(self):
        return numbers.weighed_standard_deviation_numpy(self.values, self.probabilities, mean=self.continuous_most_probable_value)

# -----------------------------------------------------------------
