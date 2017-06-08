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
from ...evolve.optimize.tables import ElitismTable, CrossoverTable
from ...evolve.optimize.stepwise import load_population
from ...core.tools.serialization import load_dict
from ...core.basics.configurable import load_input
from ...core.basics.configuration import Configuration

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

    @property
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
