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

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ...core.launch.timing import TimingTable
from ...core.launch.memory import MemoryTable
from .tables import GenerationsTable, ChiSquaredTable, ParametersTable, BestParametersTable
from ...core.basics.distribution import Distribution
from ..core.model import RTModel
from ...core.simulation.skifile import SkiFile
from ...core.simulation.skifile import SkiFile
from ...core.simulation.simulation import SkirtSimulation
from .tables import ModelProbabilitiesTable
from ...core.basics.configuration import Configuration
from ...core.filter.filter import parse_filter
from ..build.representation import Representation
from ..build.component import get_representation_path
from ...core.tools.serialization import load_dict
from .tables import IndividualsTable
from ...core.tools import types, numbers
from ...evolve.analyse.database import get_best_individual_key_and_score_all_generations
from .generation import GenerationInfo, Generation
from .evaluate import get_parameter_values_from_genome
from .platform import GenerationPlatform
from ..config.parameters import parsing_types_for_parameter_types
from ..build.component import get_representation
from ...core.tools import sequences
from ..build.component import get_model_definition
from ...core.tools.utils import lazyproperty
from ...core.simulation.wavelengthgrid import WavelengthGrid
from ...core.basics.range import QuantityRange
from ...core.tools import strings
from ..basics.properties import GalaxyProperties

# -----------------------------------------------------------------

# Files
engine_filename = "engine.pickle"
prng_filename = "prng.pickle"
optimizer_config_filename = "optimizer.cfg"
configuration_filename = "configuration.cfg"
template_filename = "template.ski"
grids_filename = "grids.dat"
weights_filename = "weights.dat"
timing_filename = "timing.dat"
memory_filename = "memory.dat"
generations_filename = "generations.dat"
best_parameters_filename = "best_parameters.dat"
input_maps_filename = "input_maps.dat"

# Directories
generations_dirname = "generations"
wavelength_grids_dirname = "wavelength grids"
best_dirname = "best"
prob_dirname = "prob"
geometries_dirname = "geometries"
refitting_dirname = "refitting"
distributions_dirname = "distributions"
parameters_dirname = "parameters"

# -----------------------------------------------------------------

class FittingRun(object):
    
    """
    This class...
    """

    def __init__(self, modeling_path, name, model_name, passive=False):

        """
        The constructor ...
        :param modeling_path:
        :param name:
        :param model_name:
        :param passive:
        :return:
        """

        # Set the name of the fitting run
        self.name = name

        # Determine the fit path
        fit_path = fs.join(modeling_path, "fit")

        # Set the path for this fitting run
        if passive: self.path = fs.join(fit_path, self.name)
        else: self.path = fs.create_directory_in(fit_path, self.name)

        # Set the name of the model used
        self.model_name = model_name

        ## Optimizer:

        # Set the path to the main genetic engine
        self.main_engine_path = fs.join(self.path, engine_filename)

        # Set the path to the main PRNG
        self.main_prng_path = fs.join(self.path, prng_filename)

        # Set the path to the optimizer configuration
        self.optimizer_config_path = fs.join(self.path, optimizer_config_filename)

        ##

        # Set the path to the fitting configuration file
        self.fitting_configuration_path = fs.join(self.path, configuration_filename)

        # Set the path to the template ski file
        self.template_ski_path = fs.join(self.path, template_filename)

        # Set the path to the fit/generations directory
        if passive: self.generations_path = fs.join(self.path, generations_dirname)
        else: self.generations_path = fs.create_directory_in(self.path, generations_dirname)

        # Set the path to the fit/wavelength grids directory
        if passive: self.wavelength_grids_path = fs.join(self.path, wavelength_grids_dirname)
        else: self.wavelength_grids_path = fs.create_directory_in(self.path, wavelength_grids_dirname)

        # Set the path to the wavelength grids table
        self.wavelength_grids_table_path = fs.join(self.wavelength_grids_path, grids_filename)

        # Set the path to the fit/best directory
        if passive: self.best_path = fs.join(self.path, best_dirname)
        else: self.best_path = fs.create_directory_in(self.path, best_dirname)

        # Set the path to the fit/prob directory
        if passive: self.prob_path = fs.join(self.path, prob_dirname)
        else: self.prob_path = fs.create_directory_in(self.path, prob_dirname)

        # Set the path to the fit/geometries directory
        if passive: self.geometries_path = fs.join(self.path, geometries_dirname)
        else: self.geometries_path = fs.create_directory_in(self.path, geometries_dirname)

        # -----------------------------------------------------------------

        ## REFITTING DIRECTORY
        if passive: self.refitting_path = fs.join(self.path, refitting_dirname)
        else: self.refitting_path = fs.create_directory_in(self.path, refitting_dirname)

        ## WEIGHTS TABLE

        # Set the path to the weights table file
        self.weights_table_path = fs.join(self.path, weights_filename)

        ## TIMING TABLE

        # Set the path to the timing table file
        self.timing_table_path = fs.join(self.path, timing_filename)

        # Initialize the timing table if necessary
        if not passive and not fs.is_file(self.timing_table_path):
            timing_table = TimingTable()
            timing_table.saveto(self.timing_table_path)

        ## MEMORY TABLE

        # Set the path to the memory table file
        self.memory_table_path = fs.join(self.path, memory_filename)

        # Initialize the memory table if necessary
        if not passive and not fs.is_file(self.memory_table_path):
            memory_table = MemoryTable()
            memory_table.saveto(self.memory_table_path)

        ## GENERATIONS TABLE

        # Set the path to the generations table
        self.generations_table_path = fs.join(self.path, generations_filename)

        # Initialize the generations table if necessary
        if not passive and not fs.is_file(self.generations_table_path) and self.free_parameter_labels is not None:
            generations_table = GenerationsTable(parameters=self.free_parameter_labels, units=self.parameter_units)
            generations_table.saveto(self.generations_table_path)

        ## PROBABILITY DISTRIBUTION TABLES

        # The directory with the probability distributions for the different generations
        if passive: self.prob_generations_path = fs.join(self.prob_path, generations_dirname)
        else: self.prob_generations_path = fs.create_directory_in(self.prob_path, generations_dirname)

        # The directory with the probability distributions for the different free parameters
        if passive: self.prob_distributions_path = fs.join(self.prob_path, distributions_dirname)
        else: self.prob_distributions_path = fs.create_directory_in(self.prob_path, distributions_dirname)

        # The directory with the combined probability tables for the different free parameters
        if passive: self.prob_parameters_path = fs.join(self.prob_path, parameters_dirname)
        else: self.prob_parameters_path = fs.create_directory_in(self.prob_path, parameters_dirname)

        ## BEST PARAMETERS TABLE

        # Set the path to the best parameters table
        self.best_parameters_table_path = fs.join(self.path, best_parameters_filename)

        # Initialize the best parameters table if necessary
        if not passive and not fs.is_file(self.best_parameters_table_path) and self.free_parameter_labels is not None:
            best_parameters_table = BestParametersTable(parameters=self.free_parameter_labels, units=self.parameter_units)
            best_parameters_table.saveto(self.best_parameters_table_path)

        ## INPUT MAP PATHS FILE
        self.input_maps_file_path = fs.join(self.path, input_maps_filename)

    # -----------------------------------------------------------------

    @lazyproperty
    def weights(self):

        """
        This function ...
        :return:
        """

        # Load the weights table
        from .tables import WeightsTable
        return WeightsTable.from_file(self.weights_table_path)

    # -----------------------------------------------------------------

    @classmethod
    def from_path(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        name = fs.name(path)
        fit_path = fs.directory_of(path)
        modeling_path = fs.directory_of(fit_path)

        # Create
        return cls.from_name(modeling_path, name)

    # -----------------------------------------------------------------

    @classmethod
    def from_name(cls, modeling_path, name):

        """
        This function ...
        :param modeling_path:
        :param name:
        :return:
        """

        from .context import FittingContext
        context = FittingContext.from_modeling_path(modeling_path)

        # Get model name
        model_name = context.get_model_for_run(name)

        # Create and return
        return cls(modeling_path, name, model_name)

    # -----------------------------------------------------------------

    @property
    def fit_path(self):
        return fs.directory_of(self.path)

    # -----------------------------------------------------------------

    @property
    def modeling_path(self):
        return fs.directory_of(self.fit_path)

    # -----------------------------------------------------------------

    @property
    def object_name(self):
        return fs.name(self.modeling_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def galaxy_properties_path(self):
        from ..core.environment import properties_name, data_name
        return fs.join(self.modeling_path, data_name, properties_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def galaxy_properties(self):
        # Load the properties
        return GalaxyProperties.from_file(self.galaxy_properties_path)

    # -----------------------------------------------------------------

    @property
    def galaxy_distance(self):
        return self.galaxy_properties.distance

    # -----------------------------------------------------------------

    @property
    def galaxy_center(self):
        return self.galaxy_properties.center

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_configuration(self):
        return Configuration.from_file(self.fitting_configuration_path) if fs.is_file(self.fitting_configuration_path) else None

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_method(self):
        return self.fitting_configuration.method

    # -----------------------------------------------------------------

    @lazyproperty
    def spectral_convolution(self):
        return self.fitting_configuration.spectral_convolution

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_filters(self):
        if self.fitting_filter_names is None: return None
        else: return map(parse_filter, self.fitting_filter_names)

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_filter_names(self):
        return self.fitting_configuration.filters if self.fitting_configuration is not None else None

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_wavelengths(self):
        return [fltr.wavelength for fltr in self.fitting_filters]

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_wavelength_range(self):

        """
        This function ...
        :return:
        """

        # Create initial wavelength range
        wavelength_range = QuantityRange.infinitesimal(self.fitting_wavelengths[0])

        # Adapt for the filters
        for fltr in self.fitting_filters: wavelength_range.adjust(fltr.wavelength)

        # Return the wavelength range
        return wavelength_range

    # -----------------------------------------------------------------

    @lazyproperty
    def absolute_fitting_wavelength_range(self):

        """
        This function ...
        :return:
        """

        # Create initial wavelength range
        wavelength_range = QuantityRange.infinitesimal(self.fitting_wavelengths[0])

        # Adapt for the filters
        for fltr in self.fitting_filters: wavelength_range.adjust(fltr.range)

        # Return the wavelength range
        return wavelength_range

    # -----------------------------------------------------------------

    @lazyproperty
    def free_parameter_labels(self):

        """
        THIS FUNCTION GUARANTEES THAT THE LABELS ARE ALWAYS ORDERED ALPHABETICALLY !!
        :return:
        """

        return sorted(self.fitting_configuration.free_parameters) if self.fitting_configuration is not None else None

    # -----------------------------------------------------------------

    def index_for_parameter(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        return self.free_parameter_labels.index(label)

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
        return self.fitting_configuration.descriptions

    # -----------------------------------------------------------------

    @lazyproperty
    def parameter_units(self):
        return self.fitting_configuration.units

    # -----------------------------------------------------------------

    @lazyproperty
    def parameter_ndigits(self):
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

    @lazyproperty
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

    @lazyproperty
    def ndigits_dict(self):
        return self.parameter_ndigits

    # -----------------------------------------------------------------

    @lazyproperty
    def nbits_list(self):

        """
        This function ...
        :return: 
        """

        # NEW: EXPERIMENTAL:
        nbits_list = []
        for index in range(len(self.ndigits_list)):
            ndigits = self.ndigits_list[index]
            low = self.parameter_minima_scalar[index]
            high = self.parameter_maxima_scalar[index]
            nbits = numbers.nbits_for_ndigits_experimental(ndigits, low, high)
            nbits_list.append(nbits)

        return nbits_list

    # -----------------------------------------------------------------

    @lazyproperty
    def nbits_dict(self):

        """
        This fucntion ...
        :return: 
        """

        keys = self.free_parameter_labels
        values = self.nbits_list
        return dict(zip(keys, values))

    # -----------------------------------------------------------------

    @lazyproperty
    def genetic_settings(self):
        return self.fitting_configuration.genetic

    # -----------------------------------------------------------------

    @lazyproperty
    def random_seed(self):
        return self.genetic_settings.seed

    # -----------------------------------------------------------------

    @lazyproperty
    def grid_settings(self):
        return self.fitting_configuration.grid

    # -----------------------------------------------------------------

    @lazyproperty
    def model_definition(self):
        return get_model_definition(self.modeling_path, self.model_name)

    # -----------------------------------------------------------------

    @property
    def normalization_filters(self):
        return self.model_definition.normalization_filters

    # -----------------------------------------------------------------

    @property
    def normalization_wavelengths(self):
        return self.model_definition.normalization_wavelengths

    # -----------------------------------------------------------------

    @property
    def normalization_center_wavelengths(self):
        return self.model_definition.normalization_center_wavelengths

    # -----------------------------------------------------------------

    @property
    def normalization_effective_wavelengths(self):
        return self.model_definition.normalization_effective_wavelengths

    # -----------------------------------------------------------------

    @property
    def normalization_pivot_wavelengths(self):
        return self.model_definition.normalization_pivot_wavelengths

    # -----------------------------------------------------------------

    @property
    def normalization_mean_wavelengths(self):
        return self.model_definition.normalization_mean_wavelengths

    # -----------------------------------------------------------------

    @property
    def normalization_peak_wavelengths(self):
        return self.model_definition.normalization_peak_wavelengths

    # -----------------------------------------------------------------

    @lazyproperty
    def initial_representation_name(self):
        return self.fitting_configuration.initial_representation

    # -----------------------------------------------------------------

    @lazyproperty
    def initial_representation(self):
        name = self.initial_representation_name
        representation_path = get_representation_path(self.modeling_path, name)
        return Representation(name, self.model_name, representation_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def parameters_config(self):
        return self.fitting_configuration.parameters

    # -----------------------------------------------------------------

    @lazyproperty
    def descriptions_config(self):
        return self.fitting_configuration.descriptions

    # -----------------------------------------------------------------

    @lazyproperty
    def types_config(self):
        return self.fitting_configuration.types

    # -----------------------------------------------------------------

    @lazyproperty
    def units_config(self):
        return self.fitting_configuration.units

    # -----------------------------------------------------------------

    @lazyproperty
    def ndigits_config(self):
        return self.fitting_configuration.ndigits

    # -----------------------------------------------------------------

    @lazyproperty
    def ranges_config(self):
        return self.fitting_configuration.ranges

    # -----------------------------------------------------------------

    @lazyproperty
    def filters_config(self):
        return self.fitting_configuration.filters

    # -----------------------------------------------------------------

    @lazyproperty
    def genetic_config(self):
        return self.fitting_configuration.genetic

    # -----------------------------------------------------------------

    @lazyproperty
    def grid_config(self):
        return self.fitting_configuration.grid

    # -----------------------------------------------------------------

    @lazyproperty
    def parameter_parsing_types(self):

        """
        This function ...
        :return:
        """

        parsing_types = dict()

        # Add the options for the ranges
        for label in self.free_parameter_labels:

            # Get the parsing type
            parsing_type = parsing_types_for_parameter_types[self.types_config[label]]

            # Add to dictinoary
            parsing_types[label] = parsing_type

        # Return
        return parsing_types

    # -----------------------------------------------------------------

    @lazyproperty
    def parameter_base_types(self):

        """
        This function ...
        :return:
        """

        base_types = dict()

        for label in self.parameter_parsing_types:

            parsing_type = self.parameter_parsing_types[label]

            if parsing_type.endswith("quantity"): base_type = "real"
            elif parsing_type == "angle": base_type = "real"
            elif parsing_type == "real": base_type = "real"
            elif parsing_type == "integer": base_type = "integer"
            else: raise ValueError("Type for parameter '" + label + "' not recognized: " + parsing_type)

            # Set type
            base_types[label] = base_type

        # Return the base types
        return base_types

    # -----------------------------------------------------------------
    # END
    # -----------------------------------------------------------------

    @property
    def needs_input(self):
        # Get the input file names
        return self.ski_template.needs_input

    # -----------------------------------------------------------------

    @property
    def wavelength_grid_names(self):
        return fs.files_in_path(self.wavelength_grids_path, extension="dat", not_contains="grids", returns="name")

    # -----------------------------------------------------------------

    @property
    def has_wavelength_grids(self):
        return len(self.wavelength_grid_names) > 0

    # -----------------------------------------------------------------

    def get_wavelength_grid_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.wavelength_grids_path, name + ".dat")

    # -----------------------------------------------------------------

    def get_wavelength_grid(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get the path
        path = self.get_wavelength_grid_path(name)

        # Load the wavelength grid and return it
        return WavelengthGrid.from_skirt_input(path)

    # -----------------------------------------------------------------

    @property
    def wavelength_grid_data_names(self):

        """
        This function ...
        :return:
        """

        return fs.directories_in_path(self.wavelength_grids_path, returns="name")

    # -----------------------------------------------------------------

    def get_wavelength_grid_data_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.wavelength_grids_path, name)

    # -----------------------------------------------------------------

    def get_wavelength_grid_plot_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get data directory path
        data_path = self.get_wavelength_grid_data_path(name)

        # Determine and return the plot path
        return fs.join(data_path, "grid.pdf")

    # -----------------------------------------------------------------

    def get_wavelength_grid_data_grid_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get data directory path
        data_path = self.get_wavelength_grid_data_path(name)

        # Determine and return the grid filepath
        return fs.join(data_path, "grid.dat")

    # -----------------------------------------------------------------

    def get_wavelength_grid_data_grid(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get the path
        path = self.get_wavelength_grid_data_grid_path(name)

        # Return
        return WavelengthGrid.from_file(path)

    # -----------------------------------------------------------------

    def get_wavelength_grid_subgrid_paths(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get data directory path
        data_path = self.get_wavelength_grid_data_path(name)

        # Get filepaths
        filepaths = fs.files_in_path(data_path, extension="dat", startswith="subgrid_")

        # Create dictionary
        paths = dict()
        for filepath in filepaths:
            subgrid = fs.strip_extension(fs.name(filepath)).split("subgrid_")[1]
            paths[subgrid] = filepath

        # Return the dictionary of paths
        return paths

    # -----------------------------------------------------------------

    def get_wavelength_grid_filter_wavelengths_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get data directory path
        data_path = self.get_wavelength_grid_data_path(name)

        # Return the path
        return fs.join(data_path, "filter_wavelengths.dat")

    # -----------------------------------------------------------------

    def get_wavelength_grid_line_wavelengths_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get data directory path
        data_path = self.get_wavelength_grid_data_path(name)

        # Return the path
        return fs.join(data_path, "line_wavelengths.dat")

    # -----------------------------------------------------------------

    def get_wavelength_grid_new_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get data directory path
        data_path = self.get_wavelength_grid_data_path(name)

        # Return the path
        return fs.join(data_path, "new.dat")

    # -----------------------------------------------------------------

    def get_wavelength_grid_replaced_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get data directory path
        data_path = self.get_wavelength_grid_data_path(name)

        # Return the path
        return fs.join(data_path, "replaced.dat")

    # -----------------------------------------------------------------

    @lazyproperty
    def timing_table(self):
        return TimingTable.from_file(self.timing_table_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def memory_table(self):
        return MemoryTable.from_file(self.memory_table_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def generations_table(self):
        return GenerationsTable.from_file(self.generations_table_path)

    # -----------------------------------------------------------------

    def index_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return self.generations_table.index_for_generation(generation_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def genetic_generation_indices(self):
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
        return self.generations_table.generation_names

    # -----------------------------------------------------------------

    def get_initial_generation_name(self):

        """
        This function ...
        :return:
        """

        return "initial"

    # -----------------------------------------------------------------

    def get_genetic_generation_name(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        if index == -1: return self.get_initial_generation_name()
        return str("Generation" + str(index))

    # -----------------------------------------------------------------

    def get_genetic_generation_index(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return int(name.split("Generation")[1])

    # -----------------------------------------------------------------

    def get_previous_genetic_generation_name(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if name == self.get_initial_generation_name(): raise ValueError("There is no generation previous to the initial generation")
        index = self.get_genetic_generation_index(name)
        if index == 0: return self.get_initial_generation_name()
        else: return self.get_genetic_generation_name(index-1)

    # -----------------------------------------------------------------

    def get_next_genetic_generation_name(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if name == self.get_initial_generation_name(): index = -1
        else: index = self.get_genetic_generation_index(name)
        return self.get_genetic_generation_name(index+1)

    # -----------------------------------------------------------------

    @lazyproperty
    def grid_generations(self):
        return self.generations_table.grid_generations

    # -----------------------------------------------------------------

    @lazyproperty
    def genetic_generations(self):
        return self.generations_table.genetic_generations

    # -----------------------------------------------------------------

    @lazyproperty
    def genetic_generations_with_initial(self):
        return self.generations_table.genetic_generations_with_initial

    # -----------------------------------------------------------------

    @property
    def finished_generations(self):
        return self.generations_table.finished_generations

    # -----------------------------------------------------------------

    @property
    def nfinished_generations(self):
        return len(self.finished_generations)

    # -----------------------------------------------------------------

    @property
    def has_finished_generations(self):
        return self.nfinished_generations > 0

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

    def get_generation_info_path(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        generation_path = self.get_generation_path(generation_name)
        return fs.join(generation_path, "info.dat")

    # -----------------------------------------------------------------

    def is_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return generation_name in self.generation_names

    # -----------------------------------------------------------------

    def get_generation_path(self, generation_name):

        """
        Thi function ...
        :param generation_name:
        :return:
        """

        return fs.join(self.generations_path, generation_name)

    # -----------------------------------------------------------------

    def get_generation_info(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return GenerationInfo.from_file(self.get_generation_info_path(generation_name))

    # -----------------------------------------------------------------

    def get_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return Generation.from_path(self.get_generation_path(generation_name))

    # -----------------------------------------------------------------

    def get_generation_platform(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return GenerationPlatform.from_path(self.get_generation_path(generation_name))

    # -----------------------------------------------------------------

    def genetic_engine_path_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        generation = self.get_generation(generation_name)
        return generation.engine_path

    # -----------------------------------------------------------------

    def prng_path_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        generation = self.get_generation(generation_name)
        return generation.prng_path

    # -----------------------------------------------------------------

    def individuals_table_path_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name: 
        :return: 
        """

        generation = self.get_generation(generation_name)
        return generation.individuals_table_path

    # -----------------------------------------------------------------

    def individuals_table_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name: 
        :return: 
        """

        generation = self.get_generation(generation_name)
        return generation.individuals_table

    # -----------------------------------------------------------------

    def chi_squared_table_path_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        generation = self.get_generation(generation_name)
        return generation.chi_squared_table_path

    # -----------------------------------------------------------------

    def chi_squared_table_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        generation = self.get_generation(generation_name)
        return generation.chi_squared_table

    # -----------------------------------------------------------------

    def parameters_table_path_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        generation = self.get_generation(generation_name)
        return generation.parameters_table_path

    # -----------------------------------------------------------------

    def parameters_table_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        generation = self.get_generation(generation_name)
        return generation.parameters_table

    # -----------------------------------------------------------------

    def elitism_table_path_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        generation = self.get_generation(generation_name)
        return generation.elitism_table_path

    # -----------------------------------------------------------------

    def elitism_table_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        generation = self.get_generation(generation_name)
        return generation.elitism_table

    # -----------------------------------------------------------------

    def recurrence_path_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        generation = self.get_generation(generation_name)
        return generation.recurrence_path

    # -----------------------------------------------------------------

    def recurrence_table_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        generation = self.get_generation(generation_name)
        return generation.recurrence_table

    # -----------------------------------------------------------------

    def crossover_table_path_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        generation = self.get_generation(generation_name)
        return generation.crossover_table_path

    # -----------------------------------------------------------------

    def crossover_table_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        generation = self.get_generation(generation_name)
        return generation.crossover_table

    # -----------------------------------------------------------------

    @lazyproperty
    def ski_template(self):

        """
        This function ...
        :return:
        """

        return SkiFile(self.template_ski_path)

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
        #print(labeled_values)
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
    def statistics(self):

        """
        Thisf unction ...
        :return:
        """

        return self.context.get_statistics_for_run(self.name)

    # -----------------------------------------------------------------

    @lazyproperty
    def populations(self):

        """
        This function ...
        :return:
        """

        #from .component import get_populations

        # Load the populations data and get only for this fitting run
        #return get_populations(self.modeling_path)[self.name]
        #return self.context.populations[self.name]
        return self.context.get_populations_for_run(self.name)

    # -----------------------------------------------------------------

    def get_initial_population(self):

        """
        This function ...
        :return:
        """

        return self.populations[0]

    # -----------------------------------------------------------------

    def get_population_for_generation(self, generation_index_or_name):

        """
        This function ...
        :param generation_index_or_name:
        :return:
        """

        if generation_index_or_name == self.get_initial_generation_name(): return self.get_initial_population()

        if types.is_string_type(generation_index_or_name): generation_index = self.get_genetic_generation_index(generation_index_or_name)
        elif types.is_integer_type(generation_index_or_name): generation_index = generation_index_or_name
        else: raise ValueError("Argument must be generation index or name")

        # Return the population for this generation
        return self.populations[generation_index + 1]

    # -----------------------------------------------------------------

    @lazyproperty
    def context(self):

        """
        Thisf unction ...
        :return:
        """

        from .context import FittingContext
        return FittingContext.from_modeling_path(self.modeling_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def best_parameter_values_and_chi_squared(self):

        """
        This function ...
        :return: 
        """

        # NEW: THIS FUNCTION WAS CREATED BECAUSE RECURRENCE WAS IMPLEMENTED: THIS MEANS THAT OUR OWN TABLES
        # (THOSE WHO CONTAIN ONLY MODELS THAT HAVE TO BE SIMULATED AND HAVE NOT OCCURED AND SCORED BEFORE)

        #from .component import get_database_path

        # Get path
        database_path = self.context.database_path

        # Get generation and individual
        generation_index, individual_key, chi_squared = get_best_individual_key_and_score_all_generations(database_path, self.name, minmax="min")
        generation_index -= 1

        # Determine generation name
        generation_name = self.get_genetic_generation_name(generation_index)

        # Load the generation
        generation = self.get_generation(generation_name)

        # Get the parameter values for the generation and individual
        individuals_generation = self.get_population_for_generation(generation_name)

        # Get genome
        genome = individuals_generation[individual_key]

        # Get parameter values from genome
        # genome, fitting_run, minima, maxima, nbits, parameter_scales, gray=False
        values = get_parameter_values_from_genome(genome, self, self.parameter_minima_scalar, self.parameter_maxima_scalar, self.nbits_list, generation.parameter_scales, self.genetic_settings.gray_code)

        # OR: USE NEW GENERATION PLATFORM CLASS

        # Return the dictionary
        return values, chi_squared

    # -----------------------------------------------------------------

    @lazyproperty
    def best_parameter_values(self):
        values, chi_squared = self.best_parameter_values_and_chi_squared
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

        # Create a RTModel object and return it
        return RTModel(self.model_definition, simulation_name=best_simulation_name, chi_squared=chi_squared,
                       free_parameter_values=parameter_values, center=self.galaxy_center)

    # -----------------------------------------------------------------

    def best_simulation_name_and_parameter_values_for_generation(self, generation_name, only_finished=True):

        """
        This function ...
        :param generation_name:
        :param only_finished:
        :return:
        """

        model = self.best_model_for_generation(generation_name, only_finished=only_finished)
        return model.simulation_name, model.parameter_values

    # -----------------------------------------------------------------

    def best_simulation_name_parameter_values_and_chi_squared_for_generation(self, generation_name, only_finished=True):

        """
        This function ...
        :param generation_name:
        :param only_finished:
        :return:
        """

        model = self.best_model_for_generation(generation_name, only_finished=only_finished)
        return model.simulation_name, model.parameter_values, model.chi_squared

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

    def get_parameter_distribution(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        # Load the probability distribution
        distribution = Distribution.from_file(self.get_parameter_distribution_path(label))

        # Normalize the distribution
        #if normalized: distribution.normalize(value=1.0, method="max")

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

    def get_probability_path_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return fs.join(self.prob_generations_path, generation_name)

    # -----------------------------------------------------------------

    def get_parameter_distribution_path_for_generation(self, label, generation_name):

        """
        This function ...
        :param label:
        :param generation_name:
        :return:
        """

        # Determine the path for the table
        return fs.join(self.get_probability_path_for_generation(generation_name), label + ".dat")

    # -----------------------------------------------------------------

    def get_parameter_distribution_for_generation(self, label, generation_name):

        """
        This function ...
        :param label:
        :param generation_name:
        :return:
        """

        # Load the distribution and return it
        return Distribution.from_file(self.get_parameter_distribution_path_for_generation(label, generation_name))

    # -----------------------------------------------------------------

    @lazyproperty
    def highest_grid_generation_index(self):

        """
        This function ...
        :return:
        """

        highest_index = -1

        # Find the highest index
        names = self.grid_generations
        for name in names:
            index_string = name.split("grid_")[1]
            if not strings.is_integer(index_string): continue
            index = int(index_string)
            if index > highest_index: highest_index = index

        # Return the highest index
        if highest_index >= 0: return highest_index
        else: return None

    # -----------------------------------------------------------------

    @lazyproperty
    def highest_grid_generation_name(self):

        """
        This function ...
        :return:
        """

        index = self.highest_grid_generation_index
        if index is None: return None
        else: return "grid_" + str(index)

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
    def last_genetic_generation_info_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.last_genetic_generation_path, "info.dat")

    # -----------------------------------------------------------------

    @lazyproperty
    def last_genetic_generation_info(self):

        """
        This function ...
        :return:
        """

        return GenerationInfo.from_file(self.last_genetic_generation_info_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def last_genetic_generation(self):

        """
        This function ...
        :return:
        """

        return Generation.from_path(self.last_genetic_generation_path)

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
    def last_genetic_or_initial_generation_info_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.last_genetic_or_initial_generation_path, "info.dat")

    # -----------------------------------------------------------------

    @lazyproperty
    def last_genetic_or_initial_generation_info(self):

        """
        This function ...
        :return:
        """

        return GenerationInfo.from_file(self.last_genetic_or_initial_generation_info_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def last_genetic_or_initial_generation(self):

        """
        This function ...
        :return:
        """

        return Generation.from_path(self.last_genetic_or_initial_generation_path)

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
    def current_spectral_convolution(self):

        """
        This function ...
        :return:
        """

        # Generations exist
        if len(self.generations_table) > 0: return self.generations_table["Spectral convolution"][-1]

        # Initial value (False)
        else: return False

    # -----------------------------------------------------------------

    @lazyproperty
    def current_use_images(self):

        """
        This function ...
        :return:
        """

        # Generations exist
        if len(self.generations_table) > 0: return self.generations_table["Use images"][-1]

        # Initial value (False)
        else: return False

    # -----------------------------------------------------------------

    @lazyproperty
    def current_wavelength_grid_name(self):

        """
        This function ...
        :return:
        """

        # Generations exist
        if len(self.generations_table) > 0: return self.generations_table["Wavelength grid name"][-1]

        # Initial value
        else: return self.lowest_basic_wavelength_grid_name

    # -----------------------------------------------------------------

    @lazyproperty
    def current_wavelength_grid_npoints(self):

        """
        Thi function ...
        :return:
        """

        return int(self.current_wavelength_grid_name.split("_")[1])

    # -----------------------------------------------------------------

    @lazyproperty
    def next_wavelength_grid_name(self):

        """
        This function ...
        :return:
        """

        # No generations yet
        if len(self.generations_table) == 0: raise ValueError("No generation yet: cannot refine")

        # Get current number of points
        current_npoints = self.current_wavelength_grid_npoints

        # Basic
        if self.is_basic_current_wavelength_grid:

            current_index = self.basic_wavelength_grid_npoints.index(current_npoints)
            next_index = current_index + 1
            next_npoints = self.basic_wavelength_grid_npoints[next_index]
            return self.get_basic_wavelength_grid_name(next_npoints)

        # Refined
        elif self.is_refined_current_wavelength_grid:

            current_index = self.refined_wavelength_grid_npoints.index(current_npoints)
            next_index = current_index + 1
            next_npoints = self.refined_wavelength_grid_npoints[next_index]
            return self.get_refined_wavelength_grid_name(next_npoints)

        # High-resolution
        elif self.is_highres_current_wavelength_grid:

            current_index = self.highres_wavelength_grid_npoints.index(current_npoints)
            next_index = current_index + 1
            next_npoints = self.highres_wavelength_grid_npoints[next_index]
            return self.get_highres_wavelength_grid_name(next_npoints)

        # Invalid
        else: raise RuntimeError("Something went wrong")

    # -----------------------------------------------------------------

    @property
    def is_basic_current_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        return self.current_wavelength_grid_name.startswith("basic")

    # -----------------------------------------------------------------

    @property
    def is_refined_current_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        return self.current_wavelength_grid_name.startswith("refined")

    # -----------------------------------------------------------------

    @property
    def is_highres_current_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        return self.current_wavelength_grid_name.startswith("highres")

    # -----------------------------------------------------------------

    @lazyproperty
    def basic_wavelength_grid_names(self):
        return fs.files_in_path(self.wavelength_grids_path, startswith="basic_", extension="dat", returns="name")

    # -----------------------------------------------------------------

    @lazyproperty
    def refined_wavelength_grid_names(self):
        return fs.files_in_path(self.wavelength_grids_path, startswith="refined_", extension="dat", returns="name")

    # -----------------------------------------------------------------

    @lazyproperty
    def highres_wavelength_grid_names(self):
        return fs.files_in_path(self.wavelength_grids_path, startswith="highres_", extension="dat", returns="name")

    # -----------------------------------------------------------------

    @lazyproperty
    def basic_wavelength_grid_npoints(self):
        return sequences.ordered([int(name.split("basic_")[1]) for name in self.basic_wavelength_grid_names])

    # -----------------------------------------------------------------

    def get_basic_wavelength_grid_name(self, npoints):

        """
        This function ...
        :param npoints:
        :return:
        """

        return "basic_" + str(npoints)

    # -----------------------------------------------------------------

    def get_basic_wavelength_grid_path(self, npoints):

        """
        This function ...
        :param npoints:
        :return:
        """

        return fs.join(self.wavelength_grids_path, self.get_basic_wavelength_grid_name(npoints) + ".dat")

    # -----------------------------------------------------------------

    def get_basic_wavelength_grid(self, npoints):

        """
        This function ...
        :param npoints:
        :return:
        """

        return WavelengthGrid.from_skirt_input(self.get_basic_wavelength_grid_path(npoints))

    # -----------------------------------------------------------------

    @property
    def lowest_basic_wavelength_grid_npoints(self):
        return min(self.basic_wavelength_grid_npoints)

    # -----------------------------------------------------------------

    @property
    def lowest_basic_wavelength_grid_name(self):
        return self.get_basic_wavelength_grid_name(self.lowest_basic_wavelength_grid_npoints)

    # -----------------------------------------------------------------

    @property
    def highest_basic_wavelength_grid_npoints(self):
        return max(self.basic_wavelength_grid_npoints)

    # -----------------------------------------------------------------

    @property
    def highest_basic_wavelength_grid_name(self):
        return self.get_basic_wavelength_grid_name(self.highest_basic_wavelength_grid_npoints)

    # -----------------------------------------------------------------

    @lazyproperty
    def refined_wavelength_grid_npoints(self):
        return sequences.ordered([int(name.split("refined_")[1]) for name in self.refined_wavelength_grid_names])

    # -----------------------------------------------------------------

    def get_refined_wavelength_grid_name(self, npoints):

        """
        This function ...
        :param npoints:
        :return:
        """

        return "refined_" + str(npoints)

    # -----------------------------------------------------------------

    def get_refined_wavelength_grid_path(self, npoints):

        """
        This function ...
        :param npoints:
        :return:
        """

        return fs.join(self.wavelength_grids_path, self.get_refined_wavelength_grid_name(npoints) + ".dat")

    # -----------------------------------------------------------------

    def get_refined_wavelength_grid(self, npoints):

        """
        This function ...
        :param npoints:
        :return:
        """

        return WavelengthGrid.from_file(self.get_refined_wavelength_grid_path(npoints))

    # -----------------------------------------------------------------

    @property
    def lowest_refined_wavelength_grid_npoints(self):
        return min(self.refined_wavelength_grid_npoints)

    # -----------------------------------------------------------------

    @property
    def lowest_refined_wavelength_grid_name(self):
        return self.get_refined_wavelength_grid_name(self.lowest_refined_wavelength_grid_npoints)

    # -----------------------------------------------------------------

    @property
    def highest_refined_wavelength_grid_npoints(self):
        return max(self.refined_wavelength_grid_npoints)

    # -----------------------------------------------------------------

    @property
    def highest_refined_wavelength_grid_name(self):
        return self.get_refined_wavelength_grid_name(self.highest_refined_wavelength_grid_npoints)

    # -----------------------------------------------------------------

    @lazyproperty
    def highres_wavelength_grid_npoints(self):
        return sequences.ordered([int(name.split("highres_")[1]) for name in self.highres_wavelength_grid_names])

    # -----------------------------------------------------------------

    def get_highres_wavelength_grid_name(self, npoints):

        """
        This function ...
        :param npoints:
        :return:
        """

        return "highres_" + str(npoints)

    # -----------------------------------------------------------------

    def get_highres_wavelength_grid_path(self, npoints):

        """
        This function ...
        :param npoints:
        :return:
        """

        return fs.join(self.wavelength_grids_path, self.get_highres_wavelength_grid_name(npoints) + ".dat")

    # -----------------------------------------------------------------

    def get_highres_wavelength_grid(self, npoints):

        """
        This function ...
        :param npoints:
        :return:
        """

        return WavelengthGrid.from_skirt_input(self.get_highres_wavelength_grid_path(npoints))

    # -----------------------------------------------------------------

    @property
    def lowest_highres_wavelength_grid_npoints(self):
        return min(self.highres_wavelength_grid_npoints)

    # -----------------------------------------------------------------

    @property
    def lowest_highres_wavelength_grid_name(self):
        return self.get_highres_wavelength_grid_name(self.lowest_highres_wavelength_grid_npoints)

    # -----------------------------------------------------------------

    @property
    def highest_highres_wavelength_grid_npoints(self):
        return max(self.highres_wavelength_grid_npoints)

    # -----------------------------------------------------------------

    @property
    def highest_highres_wavelength_grid_name(self):
        return self.get_highres_wavelength_grid_name(self.highest_highres_wavelength_grid_npoints)

    # -----------------------------------------------------------------

    @lazyproperty
    def current_model_representation_name(self):
        if len(self.generations_table) > 0: return self.generations_table["Model representation"][-1]
        else: return self.initial_representation_name #return None

    # -----------------------------------------------------------------

    @lazyproperty
    def current_model_representation(self):

        """
        This function ...
        :return:
        """

        #name = self.current_model_representation_name
        #if name is None: return self.initial_representation
        #else: return get_representation(self.modeling_path, name)

        return get_representation(self.modeling_path, self.current_model_representation_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def current_model_representation_index(self):

        """
        This function ...
        :return:
        """

        return self.model_representation_index_for_name(self.current_model_representation_name)

    # -----------------------------------------------------------------

    def model_representation_index_for_name(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return int(name.split("grid")[1])

    # -----------------------------------------------------------------

    def model_representation_name_for_index(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return "grid" + str(index)

    # -----------------------------------------------------------------

    @lazyproperty
    def previous_model_representation_index(self):

        """
        This function ...
        :return:
        """

        return self.current_model_representation_index - 1

    # -----------------------------------------------------------------

    @lazyproperty
    def previous_model_representation_name(self):

        """
        This function ...
        :return:
        """

        return self.model_representation_name_for_index(self.previous_model_representation_index)

    # -----------------------------------------------------------------

    @lazyproperty
    def previous_model_representation(self):

        """
        This function ...
        :return:
        """

        return get_representation(self.modeling_path, self.previous_model_representation_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def next_model_representation_index(self):

        """
        This function ...
        :return:
        """

        return self.current_model_representation_index + 1

    # -----------------------------------------------------------------

    @lazyproperty
    def next_model_representation_name(self):

        """
        This function ...
        :return:
        """

        return self.model_representation_name_for_index(self.next_model_representation_index)

    # -----------------------------------------------------------------

    @lazyproperty
    def next_model_representation(self):

        """
        This function ...
        :return:
        """

        return get_representation(self.modeling_path, self.next_model_representation_name)

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

class FittingRuns(object):

    """
    This function ...
    """

    def __init__(self, modeling_path):

        """
        This function ...
        :param modeling_path:
        """

        self.modeling_path = modeling_path

    # -----------------------------------------------------------------

    @lazyproperty
    def fit_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.modeling_path, "fit")

    # -----------------------------------------------------------------

    @lazyproperty
    def names(self):

        """
        This function ...
        :return:
        """

        return fs.directories_in_path(self.fit_path, returns="name")

    # -----------------------------------------------------------------

    @lazyproperty
    def empty(self):

        """
        This function ...
        :return:
        """

        return sequences.is_empty(self.names)

    # -----------------------------------------------------------------

    @lazyproperty
    def has_single(self):

        """
        This function ...
        :return:
        """

        return sequences.is_singleton(self.names)

    # -----------------------------------------------------------------

    def __len__(self):

        """
        This function ...
        :return:
        """

        return len(self.names)

    # -----------------------------------------------------------------

    @property
    def has_any(self):

        """
        This function ...
        :return:
        """

        return len(self) > 0

    # -----------------------------------------------------------------

    @lazyproperty
    def single_name(self):

        """
        This function ...
        :return:
        """

        return sequences.get_singleton(self.names)

    # -----------------------------------------------------------------

    @lazyproperty
    def single_path(self):

        """
        This function ...
        :return:
        """

        return self.get_path(self.single_name)

    # -----------------------------------------------------------------

    def get_path(self, name):

        """
        This function ...
        :param fitting_run:
        :return:
        """

        return fs.join(self.fit_path, name)

    # -----------------------------------------------------------------

    def load(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        fitting_run_path = self.get_path(name)
        if not fs.is_directory(fitting_run_path): raise ValueError("Fitting run '" + name + "' does not exist")
        return FittingRun.from_path(fitting_run_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def single(self):

        """
        This function ...
        :return:
        """

        return FittingRun.from_path(self.single_path)

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
    :param modeling_path: the modeling path
    :param fitting_run: the fitting run name
    :param generation_name: the generation name
    :return:
    """

    from .context import get_model_name_for_run

    # Get the galaxy properties
    from ..component.galaxy import get_galaxy_properties
    properties = get_galaxy_properties(modeling_path)

    # Get the model name
    model_name = get_model_name_for_run(modeling_path, fitting_run)

    # Get the model definition
    definition = get_model_definition(modeling_path, model_name)

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

    # Create a RTModel object and return it
    return RTModel(definition, simulation_name=best_simulation_name, chi_squared=chi_squared,
                   free_parameter_values=parameter_values, center=properties.center)

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

def get_fitting_configuration_path(modeling_path, fitting_run):

    """
    This function ...
    :param modeling_path: 
    :param fitting_run: 
    :return: 
    """

    #run_path = get_fitting_run_path(modeling_path, fitting_run)
    run_path = fs.join(modeling_path, "fit", fitting_run)
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
