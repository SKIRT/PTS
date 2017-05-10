#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.run Contains the FittingRun class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ...core.launch.timing import TimingTable
from ...core.launch.memory import MemoryTable
from .tables import GenerationsTable, ChiSquaredTable, ParametersTable, BestParametersTable
from ...core.simulation.skifile import LabeledSkiFile
from ...core.basics.distribution import Distribution
from ..core.model import Model
from ...core.simulation.skifile import SkiFile
from ...core.simulation.simulation import SkirtSimulation
from .tables import ModelProbabilitiesTable
from ...core.basics.configuration import Configuration
from ...core.filter.filter import parse_filter
from ..build.representation import Representation
from ..build.component import get_representation_path
from ...core.tools.serialization import load_dict
from .tables import IndividualsTable
from ...core.tools import types
from ...evolve.analyse.statistics import get_best_score_and_index_for_generation

# -----------------------------------------------------------------

class FittingRun(object):
    
    """
    This class...
    """

    def __init__(self, modeling_path, name, model_name):

        """
        The constructor ...
        :param modeling_path:
        :param name:
        :param model_name:
        :return:
        """

        # Set the name of the fitting run
        self.name = name

        # Determine the fit path
        fit_path = fs.join(modeling_path, "fit")

        # Set the path for this fitting run
        self.path = fs.create_directory_in(fit_path, self.name)

        # Set the name of the model used
        self.model_name = model_name

        ## Optimizer:

        # Set the path to the main genetic engine
        self.main_engine_path = fs.join(self.path, "engine.pickle")

        # Set the path to the main PRNG
        self.main_prng_path = fs.join(self.path, "prng.pickle")

        # Set the path to the optimizer configuration
        self.optimizer_config_path = fs.join(self.path, "optimizer.cfg")

        ##

        # Set the path to the fitting configuration file
        self.fitting_configuration_path = fs.join(self.path, "configuration.cfg")

        # Set the path to the template ski file
        self.template_ski_path = fs.join(self.path, "template.ski")

        # Set the path to the fit/generations directory
        self.generations_path = fs.create_directory_in(self.path, "generations")

        # Set the path to the fit/wavelength grids directory
        self.wavelength_grids_path = fs.create_directory_in(self.path, "wavelength grids")

        # Set the path to the wavelength grids table
        self.wavelength_grids_table_path = fs.join(self.wavelength_grids_path, "grids.dat")

        # Set the path to the fit/dust grids directory
        #self.dust_grids_path = fs.create_directory_in(self.path, "dust grids")

        # Set the path to the dust grids table
        #self.dust_grids_table_path = fs.join(self.dust_grids_path, "grids.dat")

        # Set the path to the fit/best directory
        self.best_path = fs.create_directory_in(self.path, "best")

        # Set the path to the fit/prob directory
        self.prob_path = fs.create_directory_in(self.path, "prob")

        # Set the path to the fit/instruments directory
        #self.instruments_path = fs.create_directory_in(self.path, "instruments")

        # Set the path to the SED and frame instrument
        # NOW MOVED TO REPRESENTATION
        #self.sed_instrument_path = fs.join(self.instruments_path, "sed.instr")
        #self.frame_instrument_path = fs.join(self.instruments_path, "frame.instr")
        #self.simple_instrument_path = fs.join(self.instruments_path, "simple.instr")

        # Set the path to the fit/geometries directory
        self.geometries_path = fs.create_directory_in(self.path, "geometries")

        # -----------------------------------------------------------------

        ## WEIGHTS TABLE

        # Set the path to the weights table file
        self.weights_table_path = fs.join(self.path, "weights.dat")

        ## TIMING TABLE

        # Set the path to the timing table file
        self.timing_table_path = fs.join(self.path, "timing.dat")

        # Initialize the timing table if necessary
        if not fs.is_file(self.timing_table_path):
            timing_table = TimingTable()
            timing_table.saveto(self.timing_table_path)

        ## MEMORY TABLE

        # Set the path to the memory table file
        self.memory_table_path = fs.join(self.path, "memory.dat")

        # Initialize the memory table if necessary
        if not fs.is_file(self.memory_table_path):
            memory_table = MemoryTable()
            memory_table.saveto(self.memory_table_path)

        ## GENERATIONS TABLE

        # Set the path to the generations table
        self.generations_table_path = fs.join(self.path, "generations.dat")

        # Initialize the generations table if necessary
        if not fs.is_file(self.generations_table_path) and self.free_parameter_labels is not None:
            generations_table = GenerationsTable(parameters=self.free_parameter_labels, units=self.parameter_units)
            generations_table.saveto(self.generations_table_path)

        ## PROBABILITY DISTRIBUTION TABLES

        # The directory with the probability distributions for the different free parameters
        self.prob_distributions_path = fs.create_directory_in(self.prob_path, "distributions")

        # The directory with the combined probability tables for the different free parameters
        self.prob_parameters_path = fs.create_directory_in(self.prob_path, "parameters")

        ## BEST PARAMETERS TABLE

        # Set the path to the best parameters table
        self.best_parameters_table_path = fs.join(self.path, "best_parameters.dat")

        # Initialize the best parameters table if necessary
        if not fs.is_file(self.best_parameters_table_path) and self.free_parameter_labels is not None:
            best_parameters_table = BestParametersTable(parameters=self.free_parameter_labels,
                                                        units=self.parameter_units)
            best_parameters_table.saveto(self.best_parameters_table_path)

        ## INPUT MAP PATHS FILE
        self.input_maps_file_path = fs.join(self.path, "input_maps.dat")

    # -----------------------------------------------------------------

    @property
    def modeling_path(self):

        """
        This function ...
        :return:
        """

        fit_path = fs.directory_of(self.path)
        return fs.directory_of(fit_path)

    # -----------------------------------------------------------------
    # NEW FROM MODELINGCOMPONENT
    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_configuration(self):

        """
        This function ...
        :return:
        """

        return Configuration.from_file(self.fitting_configuration_path) if fs.is_file(self.fitting_configuration_path) else None

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_method(self):

        """
        This function ...
        :return: 
        """

        return self.fitting_configuration.method

    # -----------------------------------------------------------------

    @lazyproperty
    def spectral_convolution(self):

        """
        This function ...
        :return:
        """

        return self.fitting_configuration.spectral_convolution

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_filters(self):

        """
        This function ...
        :return:
        """

        return map(parse_filter, self.fitting_filter_names)

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_filter_names(self):

        """
        This function ...
        :return:
        """

        return self.fitting_configuration.filters if self.fitting_configuration is not None else None

    # -----------------------------------------------------------------

    @lazyproperty
    def free_parameter_labels(self):

        """
        THIS FUNCTION GUARANTEES THAT THE LABELS ARE ALWAYS ORDERED ALPHABETICALLY !!
        :return:
        """

        return sorted(self.fitting_configuration.free_parameters) if self.fitting_configuration is not None else None

    # -----------------------------------------------------------------

    @lazyproperty
    def nfree_parameters(self):

        """
        This function ...
        :return:
        """

        return len(self.free_parameter_labels)

    # -----------------------------------------------------------------

    @lazyproperty
    def free_parameter_ranges(self):

        """
        This function ...
        :return:
        """

        ranges = dict()
        for label in self.free_parameter_labels:
            parameter_range = self.fitting_configuration[label + "_range"]
            ranges[label] = parameter_range
        return ranges

    # -----------------------------------------------------------------

    @lazyproperty
    def parameter_minima(self):

        """
        This function ...
        :return:
        """

        # Initialize a list
        minima = []

        # Set the list values
        for label in self.free_parameter_labels: minima.append(self.free_parameter_ranges[label].min)

        # Return the minimal parameter values
        return minima

    # -----------------------------------------------------------------

    @lazyproperty
    def parameter_maxima(self):

        """
        This function ...
        :return:
        """

        # Initialize a list
        maxima = []

        # Set the list values
        for label in self.free_parameter_labels: maxima.append(self.free_parameter_ranges[label].max)

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
        for label in self.free_parameter_labels:

            min_value = self.free_parameter_ranges[label].min

            # Convert if necessary
            if label in self.parameter_units and self.parameter_units[label] is not None:
                unit = self.parameter_units[label]
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
        for label in self.free_parameter_labels:

            max_value = self.free_parameter_ranges[label].max

            # Convert if necessary
            if label in self.parameter_units and self.parameter_units[label] is not None:
                unit = self.parameter_units[label]
                max_value = max_value.to(unit).value

            # Assert that is real type
            assert types.is_real_type(max_value)
            max_value = float(max_value)

            # Add to list
            maxima.append(max_value)

        # Return the maximal parameter values
        return maxima

    # -----------------------------------------------------------------

    @lazyproperty
    def parameter_descriptions(self):

        """
        This function ...
        :return:
        """

        return self.fitting_configuration.descriptions

    # -----------------------------------------------------------------

    @lazyproperty
    def parameter_units(self):

        """
        This function ...
        :return:
        """

        return self.fitting_configuration.units

    # -----------------------------------------------------------------

    @lazyproperty
    def parameter_ndigits(self):

        """
        This function ...
        :return: 
        """

        return self.fitting_configuration.ndigits

    # -----------------------------------------------------------------

    def ndigits_for_parameter(self, label):

        """
        This function ...
        :param label: 
        :return: 
        """

        return self.parameter_ndigits[label]

    # -----------------------------------------------------------------

    @property
    def ndigits_list(self):

        """
        This function ...
        :return: 
        """

        ndigits = []
        for label in self.free_parameter_labels:
            ndigits.append(self.ndigits_for_parameter(label))
        return ndigits

    # -----------------------------------------------------------------

    @property
    def ndigits_dict(self):

        """
        This function ...
        :return: 
        """

        return self.parameter_ndigits

    # -----------------------------------------------------------------

    @lazyproperty
    def genetic_settings(self):

        """
        This function ...
        :return:
        """

        return self.fitting_configuration.genetic

    # -----------------------------------------------------------------

    @lazyproperty
    def grid_settings(self):

        """
        This function ...
        :return: 
        """

        return self.fitting_configuration.grid

    # -----------------------------------------------------------------

    @lazyproperty
    def initial_representation(self):

        """
        This function ...
        :return:
        """

        name = self.fitting_configuration.initial_representation
        representation_path = get_representation_path(self.modeling_path, name)
        return Representation(name, self.model_name, representation_path)

    # -----------------------------------------------------------------
    # END
    # -----------------------------------------------------------------

    @property
    def needs_input(self):

        """
        This function ...
        :return:
        """

        # Get the input file names
        return self.ski_template.needs_input

    # -----------------------------------------------------------------

    @property
    def has_wavelength_grids(self):

        """
        This function ...
        :return:
        """

        return len(fs.files_in_path(self.wavelength_grids_path, extension="txt", not_contains="grids")) > 0

    # -----------------------------------------------------------------

    #@property
    #def has_dust_grids(self):

        #"""
        #This function ...
        #:return:
        #"""

        #return len(fs.files_in_path(self.dust_grids_path, extension="txt")) > 0

    # -----------------------------------------------------------------

    @lazyproperty
    def generations_table(self):

        """
        This function ...
        :return:
        """

        return GenerationsTable.from_file(self.generations_table_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def genetic_generation_indices(self):

        """
        This function ...
        :return: 
        """

        return range(self.ngenetic_generations)

    # -----------------------------------------------------------------

    @lazyproperty
    def genetic_generation_indices_for_statistics_and_database(self):

        """
        This function ...
        :return: 
        """

        return [index + 1 for index in self.genetic_generation_indices]

    # -----------------------------------------------------------------

    @lazyproperty
    def generation_names(self):

        """
        This function ...
        :return:
        """

        return self.generations_table.generation_names

    # -----------------------------------------------------------------

    @lazyproperty
    def genetic_generations(self):

        """
        This function ...
        :return:
        """

        return self.generations_table.genetic_generations

    # -----------------------------------------------------------------

    @lazyproperty
    def genetic_generations_with_initial(self):

        """
        This function ...
        :return:
        """

        return self.generations_table.genetic_generations_with_initial

    # -----------------------------------------------------------------

    @property
    def finished_generations(self):

        """
        This function ...
        :return:
        """

        return self.generations_table.finished_generations

    # -----------------------------------------------------------------

    def is_finished_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return self.generations_table.is_finished(generation_name)

    # -----------------------------------------------------------------

    @property
    def unevaluated_finished_generations(self):

        """
        This function ...
        :return:
        """

        # INitialize list
        names = []

        modeling_path = fs.directory_of(fs.directory_of(self.path))

        # Loop over the finished generations
        for generation_name in self.finished_generations:

            # Get the probabilities table
            prob_table = get_model_probabilities_table(modeling_path, self.name, generation_name)

            # If table doesn't exist yet
            if prob_table is None: names.append(generation_name)

            # Loop over all the simulation names of the generation
            for simulation_name in self.get_simulations_in_generation(generation_name):

                if not prob_table.has_simulation(simulation_name):
                    names.append(generation_name)
                    break

            else: pass # break is not encountered for this generation

        # Return the generation names
        return names

    # -----------------------------------------------------------------

    def parameter_ranges_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return self.generations_table.parameter_ranges_for_generation(generation_name)

    # -----------------------------------------------------------------

    def parameter_range_for_generation(self, generation_name, parameter_label):

        """
        This function ...
        :param generation_name: 
        :param parameter_label: 
        :return: 
        """

        return self.generations_table.parameter_range_for_generation(generation_name, parameter_label)

    # -----------------------------------------------------------------

    def parameter_minima_for_generation(self, generation_name):

        """
        THis function ...
        :param generation_name: 
        :return: 
        """

        # Initialize a list
        minima = []

        # Set the list values
        for label in self.free_parameter_labels: minima.append(self.parameter_range_for_generation(generation_name, label).min)

        # Return the minimal parameter values
        return minima

    # -----------------------------------------------------------------

    def parameter_maxima_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name: 
        :return: 
        """

        # Initialize a list
        maxima = []

        # Set the list values
        for label in self.free_parameter_labels: maxima.append(self.parameter_range_for_generation(generation_name, label).max)

        # Return the maximal parameter values
        return maxima

    # -----------------------------------------------------------------

    def parameter_minima_for_generation_scalar(self, generation_name):

        """
        This function ...
        :param generation_name: 
        :return: 
        """

        # Initialize a list
        minima = []

        # Set the list values
        for label in self.free_parameter_labels:

            # Get minimum value
            min_value = self.parameter_range_for_generation(generation_name, label).min

            # Convert if necessary
            if label in self.parameter_units and self.parameter_units[label] is not None:
                unit = self.parameter_units[label]
                min_value = min_value.to(unit).value

            # Assert that is real type
            assert types.is_real_type(min_value)
            min_value = float(min_value)

            # Add to list
            minima.append(min_value)

        # Return the minimal parameter values
        return minima

    # -----------------------------------------------------------------

    def parameter_maxima_for_generation_scalar(self, generation_name):

        """
        This function ...
        :param generation_name: 
        :return: 
        """

        # Initialize a list
        maxima = []

        # Set the list values
        for label in self.free_parameter_labels:

            max_value = self.parameter_range_for_generation(generation_name, label).max

            # Convert if necessary
            if label in self.parameter_units and self.parameter_units[label] is not None:
                unit = self.parameter_units[label]
                max_value = max_value.to(unit).value

            # Assert that is real type
            assert types.is_real_type(max_value)
            max_value = float(max_value)

            # Add to list
            maxima.append(max_value)

        # Return the maximal parameter values
        return maxima

    # -----------------------------------------------------------------

    def genetic_engine_path_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return fs.join(self.generations_path, generation_name, "engine.pickle")

    # -----------------------------------------------------------------

    def prng_path_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return fs.join(self.generations_path, generation_name, "prng.pickle")

    # -----------------------------------------------------------------

    def individuals_table_path_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name: 
        :return: 
        """

        return fs.join(self.generations_path, generation_name, "individuals.dat")

    # -----------------------------------------------------------------

    def individuals_table_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name: 
        :return: 
        """

        return IndividualsTable.from_file(self.individuals_table_path_for_generation(generation_name))

    # -----------------------------------------------------------------

    def chi_squared_table_path_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return fs.join(self.generations_path, generation_name, "chi_squared.dat")

    # -----------------------------------------------------------------

    def chi_squared_table_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return ChiSquaredTable.from_file(self.chi_squared_table_path_for_generation(generation_name))

    # -----------------------------------------------------------------

    def parameters_table_path_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return fs.join(self.generations_path, generation_name, "parameters.dat")

    # -----------------------------------------------------------------

    def parameters_table_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return ParametersTable.from_file(self.parameters_table_path_for_generation(generation_name))

    # -----------------------------------------------------------------

    @lazyproperty
    def ski_template(self):

        """
        This function ...
        :return:
        """

        return LabeledSkiFile(self.template_ski_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def first_guess_parameter_values(self):

        """
        This function ...
        :return:
        """

        # Get the current values in the ski file prepared by InputInitializer
        # young_luminosity_guess, young_filter = self.ski_template.get_stellar_component_luminosity("Young stars")
        # ionizing_luminosity_guess, ionizing_filter = self.ski_template.get_stellar_component_luminosity("Ionizing stars")
        # dust_mass_guess = self.ski_template.get_dust_component_mass(0)

        parameter_values = dict()

        # Get the values for the free parameters from the ski file template
        labeled_values = self.ski_template.get_labeled_values()
        for label in self.free_parameter_labels:
            parameter_values[label] = labeled_values[label]

        # Return the dictionary of the values for the free parameters
        return parameter_values

    # -----------------------------------------------------------------

    @property
    def has_evaluated_models(self):

        """
        This function ...
        :return:
        """

        # If there are already multiple generations, assume that there are also evaluted models for at least one
        if len(self.generation_names) > 0: return True

        # Loop over the generations, if one model has been evaluted, return True
        for generation_name in self.generation_names:

            # Check if there are 1 or more evaluted simulations in this generation
            if len(self.get_evaluated_simulations_in_generation(generation_name)) > 0: return True

        # If no evaluted simulation is found, return False
        return False

    # -----------------------------------------------------------------

    def get_evaluated_simulations_in_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        # Open the chi squared table for the specified generation
        chi_squared_table = self.chi_squared_table_for_generation(generation_name)

        # Return the names of the finished simulations
        return chi_squared_table.simulation_names

    # -----------------------------------------------------------------

    def get_simulations_in_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        # Get the parameters table
        parameters_table = self.parameters_table_for_generation(generation_name)

        # Return the names of the simulations
        return parameters_table.simulation_names

    # -----------------------------------------------------------------

    @lazyproperty
    def best_parameter_values(self):

        """
        This function ...
        :return: 
        """

        # NEW: THIS FUNCTION WAS CREATED BECAUSE RECURRENCE WAS IMPLEMENTED: THIS MEANS THAT OUR OWN TABLES
        # (THOSE WHO CONTAIN ONLY MODELS THAT HAVE TO BE SIMULATED AND HAVE NOT OCCURED AND SCORED BEFORE)

        from .component import get_statistics_path, get_populations

        # Get path
        statistics_path = get_statistics_path(self.modeling_path)

        generation_index = None
        individual_index = None
        chi_squared = float("inf")

        # Loop over the generations
        for index in self.genetic_generation_indices_for_statistics_and_database:

            # Get the best score from the statistics
            ind_index, score = get_best_score_and_index_for_generation(statistics_path, self.name, index, minmax="min")

            if score < chi_squared:
                chi_squared = score
                generation_index = index
                individual_index = ind_index

        # Look in the populations data for the parameters, for this fitting run
        populations = get_populations(self.modeling_path)[self.name]

        print("generation:", generation_index)
        print("individual:", individual_index)
        print(populations)

        # Get the parameter values for the generation and individual
        individuals_generation = populations[generation_index]
        individual_keys = individuals_generation.keys()

        print("individuals:", individual_keys)
        print("individual index:", individual_index)

        # Get parameters
        individual_key = individual_keys[individual_index]
        parameters = individuals_generation[individual_key]

        #parameters_table.parameter_values_for_simulation(best_simulation_name)

        values = dict()

        for label_index, label in enumerate(self.free_parameter_labels):

            # Add unit
            if label in self.parameter_units:
                value = parameters[label_index] * self.parameter_units[label]
            else: value = parameters[label_index]

            values[label] = value

        # Return the dictionary
        return values

    # -----------------------------------------------------------------

    @lazyproperty
    def best_parameter_values_old_implementation(self):

        """
        This function ...
        :return:
        """

        values = None
        chi_squared = float("inf")

        # Loop over the generations
        for generation_name in self.generation_names:

            generation_values, generation_chi_squared = self.best_parameter_values_for_generation(generation_name, return_chi_squared=True, only_finished=False)
            if generation_chi_squared < chi_squared:
                values = generation_values

        # Return the values
        return values

    # -----------------------------------------------------------------

    @lazyproperty
    def best_parameter_values_finished_generations(self):

        """
        This function ...
        :return:
        """

        values = None
        chi_squared = float("inf")

        # Loop over the finished generations
        for generation_name in self.finished_generations:

            generation_values, generation_chi_squared = self.best_parameter_values_for_generation(generation_name, return_chi_squared=True)
            if generation_chi_squared < chi_squared:
                values = generation_values

        # Return the values
        return values

    # -----------------------------------------------------------------

    @lazyproperty
    def best_parameters_table(self):

        """
        This function ...
        :return:
        """

        # Open the table and return it
        return BestParametersTable.from_file(self.best_parameters_table_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def best_model(self):

        """
        This function ...
        :return:
        """

        # Set best model to None initially
        best_model = None

        # Loop over the generations
        for generation_name in self.generation_names:

            # Get the best model for this generation
            model = self.best_model_for_generation(generation_name, only_finished=False)

            # Replace the best model if necessary
            if best_model is None or model.chi_squared < best_model.chi_squared: best_model = model

        # Return the best model
        return best_model

    # -----------------------------------------------------------------

    @lazyproperty
    def best_model_finished_generations(self):

        """
        This function ...
        :return:
        """

        # Set best model to None initially
        best_model = None

        # Loop over the finished generations
        for generation_name in self.finished_generations:

            # Get the best model for this generation
            model = self.best_model_for_generation(generation_name)

            # Replace the best model if necessary
            if best_model is None or model.chi_squared < best_model.chi_squared: best_model = model

        # Return the best model
        return best_model

    # -----------------------------------------------------------------

    def best_model_for_generation(self, generation_name, only_finished=True):

        """
        This function ...
        :param generation_name:
        :param only_finished:
        :return:
        """

        # Check if the generation is finished (if this is required by the caller)
        if only_finished:
            if not self.is_finished_generation(generation_name): raise RuntimeError("The generation '" + generation_name + "' is not yet finished")

        # Open the chi squared table
        chi_squared_table = self.chi_squared_table_for_generation(generation_name)

        # Get the name of the simulation with the lowest chi squared value
        best_simulation_name = chi_squared_table.best_simulation_name

        # Open the parameters table for this generation
        parameters_table = self.parameters_table_for_generation(generation_name)

        # Get the chi squared value
        chi_squared = chi_squared_table.chi_squared_for(best_simulation_name)

        # Get the parameter values
        parameter_values = parameters_table.parameter_values_for_simulation(best_simulation_name)

        # Create a 'Model' object
        model = Model()

        # Set attributes
        model.simulation_name = best_simulation_name
        model.chi_squared = chi_squared
        model.parameter_values = parameter_values

        # Return the model
        return model

    # -----------------------------------------------------------------

    def best_parameter_values_for_generation(self, generation_name, return_chi_squared=False, only_finished=True):

        """
        This function ...
        :param generation_name:
        :param return_chi_squared:
        :param only_finished:
        :return:
        """

        # Check if the generation is finished (if this is required by the caller)
        if only_finished:
            if not self.is_finished_generation(generation_name): raise RuntimeError("The generation '" + generation_name + "' is not yet finished")

        # Open the chi squared table
        chi_squared_table = self.chi_squared_table_for_generation(generation_name)

        # Get the name of the simulation with the lowest chi squared value
        best_simulation_name = chi_squared_table.best_simulation_name

        # Open the parameters table for this generation
        parameters_table = self.parameters_table_for_generation(generation_name)

        # Return the parameters of the best simulation
        if return_chi_squared:
            chi_squared = chi_squared_table.chi_squared_for(best_simulation_name)
            return parameters_table.parameter_values_for_simulation(best_simulation_name), chi_squared
        else: return parameters_table.parameter_values_for_simulation(best_simulation_name)

    # -----------------------------------------------------------------

    def get_parameter_probabilities_path(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        # Determine the path for the table
        path = fs.join(self.prob_parameters_path, label + ".dat")

        # Return the path to the table
        return path

    # -----------------------------------------------------------------

    def get_parameter_distribution_path(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        # Determine the path for the table
        path = fs.join(self.prob_distributions_path, label + ".dat")

        # Return the path to the table
        return path

    # -----------------------------------------------------------------

    def get_parameter_distribution(self, label, normalized=True):

        """
        This function ...
        :param label:
        :param normalized:
        :return:
        """

        # Load the probability distribution
        distribution = Distribution.from_file(self.get_parameter_distribution_path(label))

        # Normalize the distribution
        if normalized: distribution.normalize(value=1.0, method="max")

        # Return the distribution
        return distribution

    # -----------------------------------------------------------------

    def has_distribution(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        return fs.is_file(self.get_parameter_distribution_path(label))

    # -----------------------------------------------------------------

    @lazyproperty
    def last_genetic_generation_index(self):

        """
        This function ...
        :return:
        """

        highest_index = -1

        # Find the highest index
        for i in range(len(self.generations_table)):
            if not self.generations_table["Generation index"].mask[i]:
                index = self.generations_table["Generation index"][i]
                if index  > highest_index: highest_index = index

        # Return the highest generation index
        return highest_index

    # -----------------------------------------------------------------

    @lazyproperty
    def last_genetic_generation_name(self):

        """
        This function ...
        :return:
        """

        highest_index = -1
        name = None

        # Find the name of the generation with the highest index
        for i in range(len(self.generations_table)):
            if not self.generations_table["Generation index"].mask[i]:
                index = self.generations_table["Generation index"][i]
                if index > highest_index:
                    highest_index = index
                    name = self.generations_table["Generation name"][i]

        # Return the name of the generation with the highest index
        return name

    # -----------------------------------------------------------------

    @lazyproperty
    def last_genetic_or_initial_generation_name(self):

        """
        This function ...
        :return:
        """

        name = self.last_genetic_generation_name

        # Check whether the initial generation exists
        if name is None and "initial" in self.generations_table["Generation name"]: name = "initial"

        # Return the name
        return name

    # -----------------------------------------------------------------

    @lazyproperty
    def last_is_initial(self):

        """
        This function ...
        :return: 
        """

        return self.last_genetic_or_initial_generation_name == "initial"

    # -----------------------------------------------------------------

    @lazyproperty
    def last_genetic_generation_path(self):

        """
        This function ...
        :return:
        """

        if self.last_genetic_generation_name is None: return None
        return fs.join(self.generations_path, self.last_genetic_generation_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def last_genetic_or_initial_generation_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.generations_path, self.last_genetic_or_initial_generation_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def last_generation_name(self):

        """
        This function ...
        :return:
        """

        return self.generations_table["Generation name"][-1]

    # -----------------------------------------------------------------

    @lazyproperty
    def last_finished_generation(self):

        """
        This function ...
        :return:
        """

        # Return the name of the last finished generation
        finished_generations = self.generations_table.finished_generations
        if len(finished_generations) > 0: return finished_generations[-1]
        else: return None

    # -----------------------------------------------------------------

    @lazyproperty
    def ngenetic_generations(self):

        """
        This function ...
        :return:
        """

        return self.last_genetic_generation_index + 1

    # -----------------------------------------------------------------

    @lazyproperty
    def ngenetic_generations_with_initial(self):

        """
        This function ...
        :return: 
        """

        return self.ngenetic_generations + 1

    # -----------------------------------------------------------------

    @lazyproperty
    def ngenerations(self):

        """
        This function ...
        :return:
        """

        return len(self.generations_table)

    # -----------------------------------------------------------------

    @lazyproperty
    def current_npackages(self):

        """
        This function ...
        :return:
        """

        # Generations exist
        if len(self.generations_table) > 0: return self.generations_table["Number of photon packages"][-1]

        # Initial value
        else: return self.ski_template.packages()

    # -----------------------------------------------------------------

    @lazyproperty
    def current_selfabsorption(self):

        """
        This function ...
        :return:
        """

        # Generations exist
        if len(self.generations_table) > 0: return self.generations_table["Self-absorption"][-1]

        # Initial value
        else: return self.ski_template.dustselfabsorption()

    # -----------------------------------------------------------------

    @lazyproperty
    def current_transient_heating(self):

        """
        This function ...
        :return:
        """

        # Generations exist
        if len(self.generations_table) > 0: return self.generations_table["Transient heating"][-1]

        # Initial value
        else: return self.ski_template.transient_dust_emissivity

    # -----------------------------------------------------------------

    @lazyproperty
    def highest_wavelength_grid_level(self):

        """
        This function ...
        :return:
        """

        # Return the last filename, sorted as integers
        return int(fs.files_in_path(self.wavelength_grids_path, not_contains="grids", extension="txt", returns="name", sort=int)[-1])

    # -----------------------------------------------------------------

    @lazyproperty
    def current_wavelength_grid_level(self):

        """
        This function ...
        :return:
        """

        # Generations exist
        if len(self.generations_table) > 0: return self.generations_table["Wavelength grid level"][-1]

        # Initial value
        else: return 0

    # -----------------------------------------------------------------

    def wavelength_grid_path_for_level(self, level):

        """
        This function ...
        :param level:
        :return:
        """

        return fs.join(self.wavelength_grids_path, str(level) + ".txt")

    # -----------------------------------------------------------------

    @lazyproperty
    def current_model_representation(self):

        """
        This function ...
        :return:
        """

        if len(self.generations_table) > 0: return self.generations_table["Model representation"][-1]
        else: return None

    # -----------------------------------------------------------------

    @lazyproperty
    def input_map_paths(self):

        """
        This funtion ...
        :return:
        """

        # Load dictionary
        return load_dict(self.input_maps_file_path)

    # -----------------------------------------------------------------

    #@lazyproperty
    #def highest_dust_grid_level(self):

        #"""
        #This function ...
        #:return:
        #"""

        # Return the last filename, sorted as integers
        #return int(fs.files_in_path(self.dust_grids_path, not_contains="grids", extension="dg", returns="name", sort=int)[-1])

    # -----------------------------------------------------------------

    #@lazyproperty
    #def current_dust_grid_level(self):

        #"""
        #This function ...
        #:return:
        #"""

        # Generations exist
        #if len(self.generations_table) > 0: return self.generations_table["Dust grid level"][-1]

        # Initial value
        #else: return 0

    # -----------------------------------------------------------------

    #def dust_grid_path_for_level(self, level):

        #"""
        #This function ...
        #:param level:
        #:return:
        #"""

        #return fs.join(self.dust_grids_path, str(level) + ".dg")

    # -----------------------------------------------------------------

    #def dust_grid_for_level(self, level):

        #"""
        #This function ...
        #:param level:
        #:return:
        #"""

        # Load and return the dust grid
        #return load_grid(self.dust_grid_path_for_level(level))

    # -----------------------------------------------------------------

    #@lazyproperty
    #def sed_instrument(self):

        #"""
        #This function ...
        #:return:
        #"""

        # Check if the file exists
        #if not fs.is_file(self.sed_instrument_path): raise RuntimeError("The SED instrument file has not been created yet. Run initialize_fit first.")

        # Load the file
        #return load_instrument(self.sed_instrument_path)

    # -----------------------------------------------------------------

    #@lazyproperty
    #def frame_instrument(self):

        #"""
        #This function ...
        #:return:
        #"""

        # Check if the file exists
        #if not fs.is_file(self.frame_instrument_path): raise RuntimeError("The frame instrument file has not been created yet. Run initialize_fit first.")

        # Load the file
        #return load_instrument(self.frame_instrument_path)

    # -----------------------------------------------------------------

    #@lazyproperty
    #def simple_instrument(self):

        #"""
        #This function ...
        #:return:
        #"""

        # Check if the file exists
        #if not fs.is_file(self.simple_instrument_path): raise RuntimeError("The simple instrument file has not been created yet. Run initialize_fit first.")

        # Load the file
        #return load_instrument(self.simple_instrument_path)

# -----------------------------------------------------------------

def get_generation_names(modeling_path, fitting_run):

    """
    This function ...
    :param modeling_path:
    :param fitting_run:
    :return:
    """

    # Get the generations table
    generations_table = get_generations_table(modeling_path, fitting_run)

    # Return the generation names
    return generations_table.generation_names

# -----------------------------------------------------------------

def get_finished_generations(modeling_path, fitting_run):

    """
    This function ...
    :param modeling_path:
    :param fitting_run:
    :return:
    """

    # Get the generations table
    generations_table = get_generations_table(modeling_path, fitting_run)

    # Return the names of the finished generations
    return generations_table.finished_generations

# -----------------------------------------------------------------

def get_last_generation_name(modeling_path, fitting_run):

    """
    This function ...
    :param modeling_path:
    :param fitting_run:
    :return:
    """

    # Get the generations table
    generations_table = get_generations_table(modeling_path, fitting_run)

    # Return the name of the last generation
    if len(generations_table) > 0: return generations_table["Generation name"][-1]
    else: return None

# -----------------------------------------------------------------

def get_last_finished_generation(modeling_path, fitting_run):

    """
    This function ...
    :param modeling_path:
    :param fitting_run:
    :return:
    """

    # Get the generations table
    generations_table = get_generations_table(modeling_path, fitting_run)

    # Return the name of the last finished generation
    finished_generations = generations_table.finished_generations
    if len(finished_generations) > 0: return finished_generations[-1]
    else: return None

# -----------------------------------------------------------------

def get_generations_table(modeling_path, fitting_run):

    """
    This function ...
    :param modeling_path:
    :param fitting_run:
    :return:
    """

    # Determine the path to the generations table
    generations_table_path = fs.join(modeling_path, "fit", fitting_run, "generations.dat")

    # Load the generations table
    generations_table = GenerationsTable.from_file(generations_table_path)

    # Return the table
    return generations_table

# -----------------------------------------------------------------

def get_ngenerations(modeling_path, fitting_run):

    """
    This function ...
    :param modeling_path:
    :param fitting_run:
    :return:
    """

    # Get the table
    generations_table = get_generations_table(modeling_path, fitting_run)
    return generations_table.ngenerations

# -----------------------------------------------------------------

def get_individuals_table(modeling_path, fitting_run, generation_name):

    """
    This function ...
    :param modeling_path: 
    :param fitting_run: 
    :param generation_name: 
    :return: 
    """

    # Determine the path
    path = fs.join(modeling_path, "fit", fitting_run, "generations", generation_name, "individuals.dat")

    # Load the table
    table = IndividualsTable.from_file(path)

    # Return the table
    return table

# -----------------------------------------------------------------

def get_chi_squared_table(modeling_path, fitting_run, generation_name):

    """
    This function ...
    :param modeling_path:
    :param fitting_run:
    :param generation_name:
    :return:
    """

    # Determine the path to the chi squared table
    path = fs.join(modeling_path, "fit", fitting_run, "generations", generation_name, "chi_squared.dat")

    # Load the table
    table = ChiSquaredTable.from_file(path)

    # Return the table
    return table

# -----------------------------------------------------------------

def get_parameters_table(modeling_path, fitting_run, generation_name):

    """
    This function ...
    :param modeling_path:
    :param fitting_run:
    :param generation_name:
    :return:
    """

    # Determine the path to the parameters table
    path = fs.join(modeling_path, "fit", fitting_run, generation_name, "parameters.dat")

    # Load the table
    table = ParametersTable.from_file(path)

    # Return the table
    return table

# -----------------------------------------------------------------

def get_best_model_for_generation(modeling_path, fitting_run, generation_name):

    """
    This function ...
    :param modeling_path:
    :param fitting_run:
    :param generation_name:
    :return:
    """

    # Open the chi squared table
    chi_squared_table = get_chi_squared_table(modeling_path, fitting_run, generation_name)

    # Get the name of the simulation with the lowest chi squared value
    best_simulation_name = chi_squared_table.best_simulation_name

    # Open the parameters table for this generation
    parameters_table = get_parameters_table(generation_name, fitting_run, generation_name)

    # Get the chi squared value
    chi_squared = chi_squared_table.chi_squared_for(best_simulation_name)

    # Get the parameter values
    parameter_values = parameters_table.parameter_values_for_simulation(best_simulation_name)

    # Create a 'Model' object
    model = Model()

    # Set attributes
    model.simulation_name = best_simulation_name
    model.chi_squared = chi_squared
    model.parameter_values = parameter_values

    # Return the model
    return model

# -----------------------------------------------------------------

def get_ski_file_for_simulation(modeling_path, fitting_run, generation_name, simulation_name):

    """
    This function ...
    :param modeling_path:
    :param fitting_run:
    :param generation_name:
    :param simulation_name:
    :return:
    """

    # Get the galaxy name
    galaxy_name = fs.name(modeling_path)

    # Determine the path to the ski file
    ski_path = fs.join(modeling_path, "fit", fitting_run, "generations", generation_name, simulation_name, galaxy_name + ".ski")

    # Load and return the ski file
    return SkiFile(ski_path)

# -----------------------------------------------------------------

def get_generation_path(modeling_path, fitting_run, generation_name):

    """
    This function ...
    :param modeling_path:
    :param fitting_run:
    :param generation_name:
    :return:
    """

    generations_path = fs.join(modeling_path, "fit", fitting_run, "generations")
    return fs.join(generations_path, generation_name)

# -----------------------------------------------------------------

def get_simulation_paths(modeling_path, fitting_run, generation_name):

    """
    This function ...
    :param modeling_path:
    :param fitting_run:
    :param generation_name:
    :return:
    """

    return fs.directories_in_path(get_generation_path(modeling_path, fitting_run, generation_name))

# -----------------------------------------------------------------

def get_simulation_names(modeling_path, fitting_run, generation_name):

    """
    This function ...
    :param modeling_path:
    :param fitting_run:
    :param generation_name:
    :return:
    """

    return [fs.name(path) for path in get_simulation_paths(modeling_path, fitting_run, generation_name)]

# -----------------------------------------------------------------

def get_wavelength_grids_path(modeling_path, fitting_run):

    """
    This function ...
    :param modeling_path:
    :param fitting_run:
    :return:
    """

    return fs.join(modeling_path, "fit", fitting_run, "wavelength grids")

# -----------------------------------------------------------------

def get_simulations(modeling_path, fitting_run, generation_name):

    """
    This function ...
    :param modeling_path:
    :param generation_name:
    :param fitting_run
    :return:
    """

    # Initialize list of simulations
    simulations = []

    # Detemrine object name
    object_name = fs.name(modeling_path)

    # Loop over the simulation directories
    for simulation_path in get_simulation_paths(modeling_path, fitting_run, generation_name):

        # Get name
        simulation_name = fs.name(simulation_path)

        # Determine paths
        ski_path = fs.join(simulation_path, object_name + ".ski")
        prefix = fs.strip_extension(fs.name(ski_path))

        # Open the ski file
        ski = SkiFile(ski_path)

        # Set input file paths
        input_filenames = ski.input_files
        input_paths = []
        maps_path = fs.join(modeling_path, "maps")
        wavelength_grids_path = get_wavelength_grids_path(modeling_path, fitting_run)
        for filename in input_filenames:
            if filename.endswith(".fits"): filepath = fs.join(maps_path, filename)
            else: filepath = fs.join(wavelength_grids_path, filename)
            input_paths.append(filepath)

        # Set output path
        output_path = fs.join(simulation_path, "out")

        # Create SkirtSimulation instance
        simulation = SkirtSimulation(prefix, input_paths, output_path, ski_path, name=simulation_name)

        # Add the simulation to the list
        simulations.append(simulation)

    # Return the simulations
    return simulations

# -----------------------------------------------------------------

def has_unfinished_generations(modeling_path, fitting_run):

    """
    This function ...
    :param modeling_path:
    :param fitting_run:
    :return:
    """

    # Open the generations table
    table = get_generations_table(modeling_path, fitting_run)
    return table.has_unfinished

# -----------------------------------------------------------------

def get_model_probabilities_table(modeling_path, fitting_run, generation_name):

    """
    This function ...
    :param modeling_path:
    :param fitting_run:
    :param generation_name:
    :return:
    """

    path = fs.join(modeling_path, "fit", fitting_run, "prob", "generations", generation_name + ".dat")
    if fs.is_file(path): return ModelProbabilitiesTable.from_file(path)
    else: return None

# -----------------------------------------------------------------

def is_evaluated(modeling_path, fitting_run, generation_name):

    """
    This function ...
    :param modeling_path:
    :param fitting_run:
    :param generation_name:
    :return:
    """

    # Get the probabilities table
    prob_table = get_model_probabilities_table(modeling_path, fitting_run, generation_name)

    if prob_table is None:
        #print("prob table is None")
        return False

    # Loop over all the simulation names of the generation
    for simulation_name in get_simulation_names(modeling_path, fitting_run, generation_name):
        #print(simulation_name, prob_table.has_simulation(simulation_name))
        if not prob_table.has_simulation(simulation_name):
            #print("here")
            return False

    # No simulation encountered that was not evaluated -> OK
    return True

# -----------------------------------------------------------------

def get_evaluated_generations(modeling_path, fitting_run):

    """
    This function ...
    :param modeling_path:
    :param fitting_run:
    :return:
    """

    generation_names = []

    # Loop over the generations
    for generation_name in get_generation_names(modeling_path, fitting_run):

        if is_evaluated(modeling_path, fitting_run, generation_name): generation_names.append(generation_name)

    # Return the generation names
    return generation_names

# -----------------------------------------------------------------

def get_unevaluated_generations(modeling_path, fitting_run):

    """
    This function ...
    :param modeling_path:
    :param fitting_run:
    :return:
    """

    generation_names = []

    # Loop over the generations
    for generation_name in get_generation_names(modeling_path, fitting_run):

        if not is_evaluated(modeling_path, fitting_run, generation_name): generation_names.append(generation_name)

    # Return the generation names
    return generation_names

# -----------------------------------------------------------------

def has_unevaluated_generations(modeling_path, fitting_run):

    """
    This function ...
    :param modeling_path:
    :param fitting_run:
    :return:
    """

    # Loop over the generations
    for generation_name in get_generation_names(modeling_path, fitting_run):

        # If at least one generation is not evaluated, return False
        if not is_evaluated(modeling_path, fitting_run, generation_name): return True

    # No generation was encountered that was not completely evaluated
    return False
    
# -----------------------------------------------------------------
# NEW FROM MODELINGCOMPONENT
# -----------------------------------------------------------------

def get_fitting_run_path(modeling_path, fitting_run):

    """
    This function ...
    :param modeling_path: 
    :param fitting_run: 
    :return: 
    """

    return fs.join(modeling_path, "fit", fitting_run)

# -----------------------------------------------------------------

def get_fitting_configuration_path(modeling_path, fitting_run):

    """
    This function ...
    :param modeling_path: 
    :param fitting_run: 
    :return: 
    """

    run_path = get_fitting_run_path(modeling_path, fitting_run)
    fitting_configuration_path = fs.join(run_path, "configuration.cfg")
    return fitting_configuration_path

# -----------------------------------------------------------------

def load_fitting_configuration(modeling_path, fitting_run):

    """"
    This function ...
    :param modeling_path:
    :param fitting_run:
    :return:
    """

    # Get path
    path = get_fitting_configuration_path(modeling_path, fitting_run)

    # Determine the path to the fitting configuration file
    if not fs.is_file(path): return None

    # Open the configuration and return it
    return Configuration.from_file(path)

# -----------------------------------------------------------------

def get_free_parameter_labels(modeling_path, fitting_run):

    """
    THIS FUNCTION GUARANTEES THAT THE LABELS ARE ALWAYS ORDERED ALPHABETICALLY !!
    :return:
    """

    fitting_configuration = load_fitting_configuration(modeling_path, fitting_run)
    return sorted(fitting_configuration.free_parameters)

# -----------------------------------------------------------------

def get_parameter_descriptions(modeling_path, fitting_run):

    """
    This function ...
    :return:
    """

    fitting_configuration = load_fitting_configuration(modeling_path, fitting_run)
    return fitting_configuration.descriptions

# -----------------------------------------------------------------

def get_fitting_method(modeling_path, fitting_run):

    """
    This function ...
    :param modeling_path: 
    :param fitting_run:
    :return: 
    """

    configuration = load_fitting_configuration(modeling_path, fitting_run)
    return configuration.method

# -----------------------------------------------------------------

def get_spectral_convolution_flag(modeling_path, fitting_run):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    fitting_configuration = load_fitting_configuration(modeling_path, fitting_run)
    return fitting_configuration.spectral_convolution

# -----------------------------------------------------------------
# NEW
# -----------------------------------------------------------------
