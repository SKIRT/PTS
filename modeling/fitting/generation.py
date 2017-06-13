#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.generation Contains the Generation class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from ...core.basics.composite import SimplePropertyComposite
from ...core.tools import filesystem as fs
from ..build.component import get_model_definition, get_representation
from .tables import IndividualsTable, ParametersTable, ChiSquaredTable
from ...evolve.optimize.tables import ElitismTable, CrossoverTable, ScoresTable
from ...evolve.optimize.stepwise import load_population
from ...core.tools.serialization import load_dict
from ...core.basics.configurable import load_input
from ...core.basics.configuration import Configuration
from ...evolve.optimize.components import get_crossover, get_crossover_origins, get_genome_class, get_mutator, get_selector, get_scaling, get_initializator
from ...core.basics.range import IntegerRange, RealRange, QuantityRange
from ...core.tools import sequences

# -----------------------------------------------------------------

class GenerationInfo(SimplePropertyComposite):

    """
    This function ...
    """

    def __init__(self, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(GenerationInfo, self).__init__()

        # Define properties
        self.add_string_property("name", "name of the generation")
        self.add_integer_property("index", "index of the generation")
        self.add_string_property("method", "generation method")
        self.add_integer_property("wavelength_grid_level", "wavelength grid level")
        self.add_string_property("model_representation_name", "model representation name")
        self.add_integer_property("nsimulations", "number of simulations")
        self.add_integer_property("npackages", "number of packages")
        self.add_boolean_property("selfabsorption", "dust self-absorption enabled")
        self.add_boolean_property("transient_heating", "transient heating enabled")
        self.add_string_property("path", "generation path")

        # Set properties
        self.set_properties(kwargs)

# -----------------------------------------------------------------

class Generation(object):
    
    """
    This class...
    """
    
    def __init__(self, info):

        """
        The constructor ...
        :return:
        """

        # Set the generation info
        self.info = info

    # -----------------------------------------------------------------

    @classmethod
    def from_path(cls, directory_path):

        """
        This function ...
        :param directory_path:
        :return:
        """

        # Determine info path
        info_path = fs.join(directory_path, "info.dat")

        # Check if present
        if not fs.is_file(info_path): raise IOError("The generation info file is not present at '" + info_path + "'")

        # Load the info
        info = GenerationInfo.from_file(info_path)

        # Create the generation object and return
        return cls(info)

    # -----------------------------------------------------------------

    @property
    def info_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "info.dat")

    # -----------------------------------------------------------------

    @property
    def name(self):

        """
        This function ...
        :return:
        """

        return self.info.name

    # -----------------------------------------------------------------

    @property
    def path(self):

        """
        This function ..
        :return:
        """

        return self.info.path

    # -----------------------------------------------------------------

    @property
    def generations_path(self):

        """
        This function ...
        :return:
        """

        return fs.directory_of(self.path)

    # -----------------------------------------------------------------

    @property
    def fitting_run_path(self):

        """
        This function ...
        :return:
        """

        return fs.directory_of(self.generations_path)

    # -----------------------------------------------------------------

    @property
    def fit_path(self):

        """
        This function ...
        :return:
        """

        return fs.directory_of(self.fitting_run_path)

    # -----------------------------------------------------------------

    @property
    def modeling_path(self):

        """
        This function ...
        :return:
        """

        return fs.directory_of(self.fit_path)

    # -----------------------------------------------------------------

    @property
    def index(self):

        """
        This function ...
        :return:
        """

        return self.info.index

    # -----------------------------------------------------------------

    @property
    def method(self):

        """
        This function ...
        :return:
        """

        return self.info.method

    # -----------------------------------------------------------------

    @property
    def wavelength_grid_level(self):

        """
        This function ...
        :return:
        """

        return self.info.wavelength_grid_level

    # -----------------------------------------------------------------

    @property
    def model_representation_name(self):

        """
        This function ...
        :return:
        """

        return self.info.model_representation_name

    # -----------------------------------------------------------------

    @lazyproperty
    def model_representation(self):

        """
        This function ...
        :return:
        """

        name = self.model_representation_name
        return get_representation(self.modeling_path, name)

    # -----------------------------------------------------------------

    @property
    def model_name(self):

        """
        THis function ...
        :return:
        """

        return self.model_representation.model_name

    # -----------------------------------------------------------------

    @lazyproperty
    def model_definition(self):

        """
        THis function ...
        :return:
        """

        name = self.model_name
        return get_model_definition(self.modeling_path, name)

    # -----------------------------------------------------------------

    @property
    def nsimulations(self):

        """
        THis function ...
        :return:
        """

        return self.info.nsimulations

    # -----------------------------------------------------------------

    @property
    def npackages(self):

        """
        This function ...
        :return:
        """

        return self.info.npackages

    # -----------------------------------------------------------------

    @property
    def selfabsorption(self):

        """
        This function ...
        :return:
        """

        return self.info.selfabsorption

    # -----------------------------------------------------------------

    @property
    def transient_heating(self):

        """
        THis function ...
        :return:
        """

        return self.info.transient_heating

    # -----------------------------------------------------------------

    @property
    def individuals_table_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "individuals.dat")

    # -----------------------------------------------------------------

    @property
    def individuals_table(self):

        """
        This function ...
        :return:
        """

        return IndividualsTable.from_file(self.individuals_table_path)

    # -----------------------------------------------------------------

    @property
    def parameters_table_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "parameters.dat")

    # -----------------------------------------------------------------

    @property
    def parameters_table(self):

        """
        This function ...
        :return:
        """

        return ParametersTable.from_file(self.parameters_table_path)

    # -----------------------------------------------------------------

    @property
    def chi_squared_table_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "chi_squared.dat")

    # -----------------------------------------------------------------

    @property
    def chi_squared_table(self):

        """
        This function ...
        :return:
        """

        return ChiSquaredTable.from_file(self.chi_squared_table_path)

    # -----------------------------------------------------------------

    @property
    def newborns_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "newborns.dat")

    # -----------------------------------------------------------------

    @property
    def newborns(self):

        """
        This function ...
        :return:
        """

        if fs.is_file(self.newborns_path): return load_population(self.newborns_path)
        else: return None

    # -----------------------------------------------------------------

    @property
    def parents_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "parents.dat")

    # -----------------------------------------------------------------

    @property
    def parents(self):

        """
        This function ...
        :return:
        """

        return load_population(self.parents_path)

    # -----------------------------------------------------------------

    @property
    def scores_table_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "scores")

    # -----------------------------------------------------------------

    @property
    def scores_table(self):

        """
        This function ...
        :return:
        """

        return ScoresTable.from_file(self.scores_table_path)

    # -----------------------------------------------------------------

    @property
    def elitism_table_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "elitism.dat")

    # -----------------------------------------------------------------

    @property
    def elitism_table(self):

        """
        This function ...
        :return:
        """

        return ElitismTable.from_file(self.elitism_table_path)

    # -----------------------------------------------------------------

    @property
    def recurrent_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "recurrent.dat")

    # -----------------------------------------------------------------

    @property
    def recurrent_data(self):

        """
        This function ...
        :return:
        """

        return load_dict(self.recurrent_path)

    # -----------------------------------------------------------------

    @property
    def crossover_table_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "crossover.dat")

    # -----------------------------------------------------------------

    @lazyproperty
    def crossover_table(self):

        """
        This function ...
        :return:
        """

        return CrossoverTable.from_file(self.crossover_table_path)

    # -----------------------------------------------------------------

    @property
    def optimizer_input_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.path, "input")

    # -----------------------------------------------------------------

    @property
    def optimizer_input(self):

        """
        This function ...
        :return:
        """

        return load_input(self.optimizer_input_path)

    # -----------------------------------------------------------------

    @property
    def parameter_minima_scalar(self):

        """
        This function ...
        :return:
        """

        return self.optimizer_input["minima"]

    # -----------------------------------------------------------------

    @property
    def parameter_maxima_scalar(self):

        """
        This function ...
        :return:
        """

        return self.optimizer_input["maxima"]

    # -----------------------------------------------------------------

    def unit_for_parameter(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        return self.parameters_table.unit_for(label)

    # -----------------------------------------------------------------

    @property
    def parameter_labels(self):

        """
        This function ...
        :return:
        """

        ## SORTED JUST AS IN FITTINGRUN !!
        return sorted(self.parameters_table.parameter_labels)

    # -----------------------------------------------------------------

    @property
    def parameter_units(self):

        """
        This function ...
        :return:
        """

        units = dict()

        for label in self.parameter_labels:

            unit = self.unit_for_parameter(label)
            units[label] = unit

        # Return the units
        return units

    # -----------------------------------------------------------------

    @property
    def parameter_minima(self):

        """
        This function ...
        :return:
        """

        minima = dict()

        for minimum, label in zip(self.parameter_minima_scalar, self.parameter_labels):

            unit = self.unit_for_parameter(label)
            if unit is not None: minimum *= unit

            minima[label] = minimum

        # Return
        return minima

    # -----------------------------------------------------------------

    @property
    def parameter_maxima(self):

        """
        This function ...
        :return:
        """

        maxima = dict()

        for maximum, label in zip(self.parameter_maxima_scalar, self.parameter_labels):

            unit = self.unit_for_parameter(label)
            if unit is not None: maximum *= unit

            maxima[label] = maximum

        # Return
        return maxima

    # -----------------------------------------------------------------

    @property
    def parameter_ranges(self):

        """
        This function ...
        :return:
        """

        ranges = dict()

        minima = self.parameter_minima
        maxima = self.parameter_maxima

        for label in self.parameter_labels:

            # Get min and max value
            min_value = minima[label]
            max_value = maxima[label]

            unit = self.unit_for_parameter(label)

            # Quantity
            if unit is not None: range = QuantityRange(min_value, max_value)

            # Scalar value
            else:

                base_type = self.parameter_base_types
                if base_type == "integer": range = IntegerRange(min_value, max_value)
                elif base_type == "real": range = RealRange(min_value, max_value)
                else: raise ValueError("Unrecognized base type: " + base_type)

            # Add the range
            ranges[label] = range

        # Return
        return ranges

    # -----------------------------------------------------------------

    @property
    def parameter_scales(self):

        """
        This function ...
        :return:
        """

        return self.optimizer_input["scales"]

    # -----------------------------------------------------------------

    @property
    def parameter_ndigits(self):

        """
        This function ...
        :return:
        """

        return self.optimizer_input["ndigits"]

    # -----------------------------------------------------------------

    @property
    def ndigits_list(self):

        """
        This function ...
        :return:
        """

        return self.parameter_ndigits

    # -----------------------------------------------------------------

    @property
    def parameter_nbits(self):

        """
        This function ...
        :return:
        """

        return self.optimizer_input["nbits"]

    # -----------------------------------------------------------------

    @property
    def nbits_list(self):

        """
        This function ...
        :return:
        """

        return self.parameter_nbits

    # -----------------------------------------------------------------

    @property
    def total_nbits(self):

        """
        This function ...
        :return:
        """

        return sum(self.nbits_list)

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
    def genome_size(self):

        """
        This function ...
        :return:
        """

        # 1D
        if self.is_1d_genome:

            # 1D
            if self.list_genome: return self.nparameters
            elif self.binary_string_genome: return self.total_nbits
            else: raise ValueError("Invalid genome type")

        # 2D
        #elif self.is_2d_genome: raise ValueError("We should not get here")
        raise ValueError("We should not get here")

    # -----------------------------------------------------------------

    @property
    def elitism(self):

        """
        This function ...
        :return:
        """

        return self.optimizer_config.elitism

    # -----------------------------------------------------------------

    @property
    def check_recurrence(self):

        """
        This function ...
        :return:
        """

        return self.optimizer_config.check_recurrence

    # -----------------------------------------------------------------

    @property
    def gray_code(self):

        """
        This function ...
        :return:
        """

        return self.optimizer_config.gray_code

    # -----------------------------------------------------------------

    @property
    def genome_type(self):

        """
        This function ...
        :return:
        """

        return self.optimizer_config.genome_type

    # -----------------------------------------------------------------

    @property
    def list_genome(self):

        """
        This function ...
        :return:
        """

        return self.genome_type == "list"

    # -----------------------------------------------------------------

    @property
    def binary_string_genome(self):

        """
        This function ...
        :return:
        """

        return self.genome_type == "binary_string"

    # -----------------------------------------------------------------

    @property
    def genome_class(self):

        """
        This function ...
        :return:
        """

        #return genomes[self.genome_dimension][self.genome_type]
        return get_genome_class(self.genome_dimension, self.genome_type)

    # -----------------------------------------------------------------

    ## FIXED FOR MODELLING
    @property
    def genome_dimension(self):

        """
        This function ...
        :return:
        """

        return 1

    # -----------------------------------------------------------------

    ## FIXED FOR MODELLING
    @property
    def is_1d_genome(self):

        """
        This fucntion ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    ## FIXED FOR MODELLING
    @property
    def is_2d_genome(self):

        """
        This function ...
        :return:
        """

        return False

    # -----------------------------------------------------------------

    ## FIXED FOR MODELLING
    @property
    def min_or_max(self):

        """
        This function ...
        :return:
        """

        return "minimize"

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_run(self):

        """
        This function ...
        :return:
        """

        from .run import FittingRun
        return FittingRun.from_path(self.fitting_run_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def parameter_base_types(self):

        """
        This function ...
        :return:
        """

        types = dict()

        ranges = self.parameter_ranges

        for label in self.parameter_labels:

            the_range = ranges[label]

            # Check the type of range
            if isinstance(the_range, IntegerRange): base_type = "integer"
            elif isinstance(the_range, RealRange): base_type = "real"
            elif isinstance(the_range, QuantityRange): base_type = "real" # quantity -> 'real' base type
            else: raise ValueError("Invalid parameter range")

            # Set the type
            types[label] = base_type

        # Return
        return types

    # -----------------------------------------------------------------

    @lazyproperty
    def single_parameter_base_type(self):

        """
        This function ...
        :return:
        """

        assert sequences.all_equal(self.parameter_base_types.values())
        return self.parameter_base_types[self.parameter_base_types.keys()[0]]

    # -----------------------------------------------------------------

    ## FIXED FOR MODELLING
    @property
    def heterogeneous(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def nindividuals_per_generation(self):

        """
        This function ...
        :return:
        """

        return self.optimizer_config.nindividuals

    # -----------------------------------------------------------------

    @property
    def mutation_rate(self):

        """
        This function ...
        :return:
        """

        return self.optimizer_config.mutation_rate

    # -----------------------------------------------------------------

    @property
    def nparameters(self):

        """
        This function ...
        :return:
        """

        return self.optimizer_config.nparameters

    # -----------------------------------------------------------------

    @property
    def crossover_rate(self):

        """
        This function ...
        :return:
        """

        return self.optimizer_config.crossover_rate

    # -----------------------------------------------------------------

    @property
    def mutation_method(self):

        """
        This function ...
        :return:
        """

        return self.optimizer_config.mutation_method

    # -----------------------------------------------------------------

    @property
    def binary_mutation_method(self):

        """
        THis function ...
        :return:
        """

        return self.optimizer_config.binary_mutation_method

    # -----------------------------------------------------------------

    @property
    def mutation_function(self):

        """
        This function ...
        :return:
        """

        # List or binary
        if self.list_genome: return get_mutator(self.genome_type, self.genome_dimension, self.mutation_method, parameter_type=self.single_parameter_base_type, hetero=self.heterogeneous)
        elif self.binary_string_genome: return get_mutator(self.genome_type, self.genome_dimension, self.binary_mutation_method)
        else: raise ValueError("Invalid genome type")

    # -----------------------------------------------------------------

    @property
    def scaling_method(self):

        """
        This function ...
        :return:
        """

        return self.optimizer_config.scaling_method

    # -----------------------------------------------------------------

    @property
    def scaling_function(self):

        """
        This function ...
        :return:
        """

        return get_scaling(self.scaling_method)

    # -----------------------------------------------------------------

    @property
    def selection_method(self):

        """
        This function ...
        :return:
        """

        return self.optimizer_config.selector_method

    # -----------------------------------------------------------------

    @property
    def selection_function(self):

        """
        This function ...
        :return:
        """

        return get_selector(self.selection_method)

    # -----------------------------------------------------------------

    @property
    def initializator(self):

        """
        This function ...
        :return:
        """

        return get_initializator(self.genome_type, parameter_type=self.single_parameter_base_type, hetero=self.heterogeneous)

    # -----------------------------------------------------------------

    @property
    def crossover_method(self):

        """
        This function ...
        :return:
        """

        return self.optimizer_config.crossover_method

    # -----------------------------------------------------------------

    @property
    def crossover_function(self):

        """
        This function ...
        :return:
        """

        return get_crossover(self.genome_type, self.genome_dimension, self.crossover_method, genome_size=self.genome_size)

    # -----------------------------------------------------------------

    @property
    def crossover_origins_function(self):

        """
        This function ...
        :return:
        """

        return get_crossover_origins(self.genome_type, self.genome_dimension, self.crossover_method)

    # -----------------------------------------------------------------

    @property
    def engine_path(self):

        """
        THis function ...
        :return:
        """

        return fs.join(self.path, "engine.pickle")

    # -----------------------------------------------------------------

    @property
    def prng_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "prng.pickle")

    # -----------------------------------------------------------------

    @property
    def optimizer_config_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "optimizer.cfg")

    # -----------------------------------------------------------------

    @property
    def optimizer_config(self):

        """
        This function ...
        :return:
        """

        return Configuration.from_file(self.optimizer_config_path)

# -----------------------------------------------------------------
