#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.component Contains the FittingComponent class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta

# Import astronomical modules
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
#from ..component.galaxy import GalaxyModelingComponent
from ..component.component import ModelingComponent
from ...core.tools import filesystem as fs
from ...core.launch.timing import TimingTable
from ...core.launch.memory import MemoryTable
from .tables import GenerationsTable, ChiSquaredTable, ParametersTable, BestParametersTable
from ...core.simulation.skifile import LabeledSkiFile
from ...core.basics.distribution import Distribution
from ..basics.instruments import load_instrument
from ..core.model import Model
from ...core.simulation.grids import load_grid
from ...core.simulation.skifile import SkiFile

# -----------------------------------------------------------------

class FittingComponent(ModelingComponent):
    
    """
    This class...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(FittingComponent, self).__init__(config)

        # -- Attributes --

        # The path to the template ski file
        self.template_ski_path = None

        # The path to the fixed parameter values file
        self.fixed_parameters_path = None

        # The path to the fit/evolution directory
        self.fit_generations_path = None

        # The path to the fit/wavelength grids directory
        self.fit_wavelength_grids_path = None

        # The path to the wavelength grids table
        self.wavelength_grids_table_path = None

        # The path to the fit/dust grids directory
        self.fit_dust_grids_path = None

        # The path to the dust grids table
        self.dust_grids_table_path = None

        # The path to the fit/best directory
        self.fit_best_path = None

        # The path to the fit/instruments directory
        self.fit_instruments_path = None

        # The paths to the SED and frame instrument files
        self.sed_instrument_path = None
        self.frame_instrument_path = None
        self.simple_instrument_path = None

        # The path to the fit/prob directory
        self.fit_prob_path = None

        # The path to the generations table
        self.generations_table_path = None

        # The path to the best parameters table
        self.best_parameters_table_path = None

        # The path to the weights table
        self.weights_table_path = None

        # The path to the timing table
        self.timing_table_path = None

        # The path to the memory table
        self.memory_table_path = None

        # The path to the geometries directory
        self.fit_geometries_path = None

        # The directory with the probability distributions for the different free parameters
        self.prob_distributions_path = None

        # The paths to the probability distribution tables
        #self.distribution_table_paths = dict()

        # The directory for the posterior / parameter probability distribution tables
        self.prob_parameters_path = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(FittingComponent, self).setup(**kwargs)

        # Set the path to the template ski file
        self.template_ski_path = fs.join(self.fit_path, "template.ski")

        # Set the path to the fixed parameters file
        self.fixed_parameters_path = fs.join(self.fit_path, "fixed.dat")

        # Set the path to the fit/generations directory
        self.fit_generations_path = fs.create_directory_in(self.fit_path, "generations")

        # Set the path to the fit/wavelength grids directory
        self.fit_wavelength_grids_path = fs.create_directory_in(self.fit_path, "wavelength grids")

        # Set the path to the wavelength grids table
        self.wavelength_grids_table_path = fs.join(self.fit_wavelength_grids_path, "grids.dat")

        # Set the path to the fit/dust grids directory
        self.fit_dust_grids_path = fs.create_directory_in(self.fit_path, "dust grids")

        # Set the path to the dust grids table
        self.dust_grids_table_path = fs.join(self.fit_dust_grids_path, "grids.dat")

        # Set the path to the fit/best directory
        self.fit_best_path = fs.create_directory_in(self.fit_path, "best")

        # Set the path to the fit/prob directory
        self.fit_prob_path = fs.create_directory_in(self.fit_path, "prob")

        # Set the path to the fit/instruments directory
        self.fit_instruments_path = fs.create_directory_in(self.fit_path, "instruments")

        # Set the path to the SED and frame instrument
        self.sed_instrument_path = fs.join(self.fit_instruments_path, "sed.instr")
        self.frame_instrument_path = fs.join(self.fit_instruments_path, "frame.instr")
        self.simple_instrument_path = fs.join(self.fit_instruments_path, "simple.instr")

        # Set the path to the fit/geometries directory
        self.fit_geometries_path = fs.create_directory_in(self.fit_path, "geometries")

        # -----------------------------------------------------------------

        ## WEIGHTS TABLE

        # Set the path to the weights table file
        self.weights_table_path = fs.join(self.fit_path, "weights.dat")

        ## TIMING TABLE

        # Set the path to the timing table file
        self.timing_table_path = fs.join(self.fit_path, "timing.dat")

        # Initialize the timing table if necessary
        if not fs.is_file(self.timing_table_path):
            timing_table = TimingTable()
            timing_table.saveto(self.timing_table_path)

        ## MEMORY TABLE

        # Set the path to the memory table file
        self.memory_table_path = fs.join(self.fit_path, "memory.dat")

        # Initialize the memory table if necessary
        if not fs.is_file(self.memory_table_path):
            memory_table = MemoryTable()
            memory_table.saveto(self.memory_table_path)

        ## GENERATIONS TABLE

        # Set the path to the generations table
        self.generations_table_path = fs.join(self.fit_path, "generations.dat")

        # Initialize the generations table if necessary
        if not fs.is_file(self.generations_table_path) and self.free_parameter_labels is not None:
            generations_table = GenerationsTable(parameters=self.free_parameter_labels, units=self.parameter_units)
            generations_table.saveto(self.generations_table_path)

        ## PROBABILITY DISTRIBUTION TABLES

        # The directory with the probability distributions for the different free parameters
        self.prob_distributions_path = fs.create_directory_in(self.fit_prob_path, "distributions")

        # The directory with the combined probability tables for the different free parameters
        self.prob_parameters_path = fs.create_directory_in(self.fit_prob_path, "parameters")

        ## BEST PARAMETERS TABLE

        # Set the path to the best parameters table
        self.best_parameters_table_path = fs.join(self.fit_path, "best_parameters.dat")

        # Initialize the best parameters table if necessary
        if not fs.is_file(self.best_parameters_table_path) and self.free_parameter_labels is not None:
            best_parameters_table = BestParametersTable(parameters=self.free_parameter_labels, units=self.parameter_units)
            best_parameters_table.saveto(self.best_parameters_table_path)

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

        return len(fs.files_in_path(self.fit_wavelength_grids_path, extension="txt")) > 0

    # -----------------------------------------------------------------

    @property
    def has_dust_grids(self):

        """
        This function ...
        :return:
        """

        return len(fs.files_in_path(self.fit_dust_grids_path, extension="txt")) > 0

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

    @lazyproperty
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

    def parameter_ranges_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return self.generations_table.parameter_ranges_for_generation(generation_name)

    # -----------------------------------------------------------------

    def genetic_engine_path_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return fs.join(self.fit_generations_path, generation_name, "engine.pickle")

    # -----------------------------------------------------------------

    def prng_path_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return fs.join(self.fit_generations_path, generation_name, "prng.pickle")

    # -----------------------------------------------------------------

    def chi_squared_table_path_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return fs.join(self.fit_generations_path, generation_name, "chi_squared.dat")

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

        return fs.join(self.fit_generations_path, generation_name, "parameters.dat")

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

    @lazyproperty
    def fixed_parameters(self):

        """
        This function ...
        :return:
        """

        from ...core.basics.configuration import load_mapping
        from ...core.basics.map import Map

        parameters = Map()
        with open(self.fixed_parameters_path) as f: load_mapping(f, parameters)

        # Return the parameters map
        return parameters

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
    def last_genetic_generation_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.fit_generations_path, self.last_genetic_generation_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def last_genetic_or_initial_generation_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.fit_generations_path, self.last_genetic_or_initial_generation_name)

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
    def ngenetic_generations(self):

        """
        This function ...
        :return:
        """

        return self.last_genetic_generation_index + 1

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
        return int(fs.files_in_path(self.fit_wavelength_grids_path, not_contains="grids", extension="txt", returns="name", sort=int)[-1])

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

        return fs.join(self.fit_wavelength_grids_path, str(level) + ".txt")

    # -----------------------------------------------------------------

    @lazyproperty
    def highest_dust_grid_level(self):

        """
        This function ...
        :return:
        """

        # Return the last filename, sorted as integers
        return int(fs.files_in_path(self.fit_dust_grids_path, not_contains="grids", extension="dg", returns="name", sort=int)[-1])

    # -----------------------------------------------------------------

    @lazyproperty
    def current_dust_grid_level(self):

        """
        This function ...
        :return:
        """

        # Generations exist
        if len(self.generations_table) > 0: return self.generations_table["Dust grid level"][-1]

        # Initial value
        else: return 0

    # -----------------------------------------------------------------

    def dust_grid_path_for_level(self, level):

        """
        This function ...
        :param level:
        :return:
        """

        return fs.join(self.fit_dust_grids_path, str(level) + ".dg")

    # -----------------------------------------------------------------

    def dust_grid_for_level(self, level):

        """
        This function ...
        :param level:
        :return:
        """

        # Load and return the dust grid
        return load_grid(self.dust_grid_path_for_level(level))

    # -----------------------------------------------------------------

    @lazyproperty
    def sed_instrument(self):

        """
        This function ...
        :return:
        """

        # Check if the file exists
        if not fs.is_file(self.sed_instrument_path): raise RuntimeError("The SED instrument file has not been created yet. Run initialize_fit first.")

        # Load the file
        return load_instrument(self.sed_instrument_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def frame_instrument(self):

        """
        This function ...
        :return:
        """

        # Check if the file exists
        if not fs.is_file(self.frame_instrument_path): raise RuntimeError("The frame instrument file has not been created yet. Run initialize_fit first.")

        # Load the file
        return load_instrument(self.frame_instrument_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def simple_instrument(self):

        """
        This function ...
        :return:
        """

        # Check if the file exists
        if not fs.is_file(self.simple_instrument_path): raise RuntimeError("The simple instrument file has not been created yet. Run initialize_fit first.")

        # Load the file
        return load_instrument(self.simple_instrument_path)

# -----------------------------------------------------------------

def get_generation_names(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    # Get the generations table
    generations_table = get_generations_table(modeling_path)

    # Return the generation names
    return generations_table.generation_names

# -----------------------------------------------------------------

def get_finished_generations(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    # Get the generations table
    generations_table = get_generations_table(modeling_path)

    # Return the names of the finished generations
    return generations_table.finished_generations

# -----------------------------------------------------------------

def get_last_generation_name(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    # Get the generations table
    generations_table = get_generations_table(modeling_path)

    # Return the name of the last generation
    if len(generations_table) > 0: return generations_table["Generation name"][-1]
    else: return None

# -----------------------------------------------------------------

def get_last_finished_generation(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    # Get the generations table
    generations_table = get_generations_table(modeling_path)

    # Return the name of the last finished generation
    finished_generations = generations_table.finished_generations
    if len(finished_generations) > 0: return finished_generations[-1]
    else: return None

# -----------------------------------------------------------------

def get_generations_table(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    # Determine the path to the generations table
    generations_table_path = fs.join(modeling_path, "fit", "generations.dat")

    # Load the generations table
    generations_table = GenerationsTable.from_file(generations_table_path)

    # Return the table
    return generations_table

# -----------------------------------------------------------------

def get_chi_squared_table(modeling_path, generation_name):

    """
    This function ...
    :param modeling_path:
    :param generation_name:
    :return:
    """

    # Determine the path to the chi squared table
    path = fs.join(modeling_path, "fit", generation_name, "chi_squared.dat")

    # Load the table
    table = ChiSquaredTable.from_file(path)

    # Return the table
    return table

# -----------------------------------------------------------------

def get_parameters_table(modeling_path, generation_name):

    """
    This function ...
    :param modeling_path:
    :param generation_name:
    :return:
    """

    # Determine the path to the parameters table
    path = fs.join(modeling_path, "fit", generation_name, )

    # Load the table
    table = ParametersTable.from_file(path)

    # Return the table
    return table

# -----------------------------------------------------------------

def get_best_model_for_generation(modeling_path, generation_name):

    """
    This function ...
    :param modeling_path:
    :param generation_name:
    :return:
    """

    # Open the chi squared table
    chi_squared_table = get_chi_squared_table(modeling_path, generation_name)

    # Get the name of the simulation with the lowest chi squared value
    best_simulation_name = chi_squared_table.best_simulation_name

    # Open the parameters table for this generation
    parameters_table = get_parameters_table(generation_name, generation_name)

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

def get_ski_file_for_simulation(modeling_path, generation_name, simulation_name):

    """
    This function ...
    :param modeling_path:
    :param generation_name:
    :param simulation_name:
    :return:
    """

    # Get the galaxy name
    galaxy_name = fs.name(modeling_path)

    # Determine the path to the ski file
    ski_path = fs.join(modeling_path, generation_name, simulation_name, galaxy_name + ".ski")

    # Load and return the ski file
    return SkiFile(ski_path)

# -----------------------------------------------------------------
