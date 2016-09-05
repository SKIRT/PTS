#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.modelgenerators.genetic Contains the GeneticModelGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....core.tools.logging import log
from ....evolve.engine import GAEngine
from .generator import ModelGenerator
from ....core.tools import filesystem as fs
from ....core.tools.random import save_state, load_state

# -----------------------------------------------------------------

class GeneticModelGenerator(ModelGenerator):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(GeneticModelGenerator, self).__init__(config)

        # -- Attributes --

        # The genetic algorithm engine
        self.engine = None

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(GeneticModelGenerator, self).setup()

        # Load the state of the prng
        self.load_random_state()

        # Load the genetic engine
        self.load_engine()

        # Set options for the genetic engine
        self.set_options()

    # -----------------------------------------------------------------

    def load_random_state(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the state of the random number generator ...")

        # Determine the path to the saved prng state
        path = fs.join(self.last_genetic_or_initial_generation_path, "prng.pickle")

        # Load the random state
        load_state(path)

    # -----------------------------------------------------------------

    def load_engine(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the state of the genetic algorithm engine ...")

        # Determine the path to the saved genetic algorithm engine
        path = fs.join(self.last_genetic_or_initial_generation_path, "engine.pickle")

        # Load the engine
        self.engine = GAEngine.from_file(path)

    # -----------------------------------------------------------------

    def set_options(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting options for the genetic algorithm engine ...")

        # Set genome ranges
        self.engine.set_ranges(minima=self.parameter_minima, maxima=self.parameter_maxima)

        # Set other options
        self.engine.setCrossoverRate(self.config.crossover_rate)
        self.engine.setPopulationSize(self.config.nmodels)
        self.engine.setMutationRate(self.config.mutation_rate)

        # Make sure we don't stop unexpectedly, always increment the number of generations we want (it doesn't matter what the value is)
        self.engine.setGenerations(self.engine.getGenerations() + 1)

    # -----------------------------------------------------------------

    def generate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the new models ...")

        # Set the scores from the previous generation
        self.set_scores()

        # Generate the new models
        self.generate_new_models()

    # -----------------------------------------------------------------

    def set_scores(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting scores from previous generation ...")

        # Load the parameters table from the previous generation
        parameters_table = self.parameters_table_for_generation(self.last_genetic_or_initial_generation_name)

        # Load the chi squared table from the previous generation
        chi_squared_table = self.chi_squared_table_for_generation(self.last_genetic_or_initial_generation_name)

        # List of chi squared values in the same order as the parameters table
        chi_squared_values = []

        # Check whether the chi-squared and parameter tables match
        for i in range(len(parameters_table)):

            simulation_name = parameters_table["Simulation name"][i]
            chi_squared = chi_squared_table.chi_squared_for(simulation_name)
            chi_squared_values.append(chi_squared)

        # Get the scores
        scores = chi_squared_table["Chi squared"]

        # Check individual values with parameter table of the last generation
        check = []
        for label in self.free_parameter_labels:
            values = parameters_table[label]
            check.append(values)

        # Set the scores
        self.engine.set_scores(scores, check)

    # -----------------------------------------------------------------

    def generate_new_models(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the new population ...")

        # Generate the new population
        self.engine.generate_new_population()

        #for ind in self.engine.new_population:

            # Give the individual a unique name
            #name = time.unique_name(precision="micro")
            #name_column.append(name)
            #par_a_column.append(ind.genomeList[0])
            #par_b_column.append(ind.genomeList[1])
            #par_c_column.append(ind.genomeList[2])

        # Loop over the individuals of the population
        for individual in self.engine.new_population:

            # Loop over all the genes (parameters)
            for i in range(len(individual)):

                # Get the parameter value
                value = individual[i]

                # Add the parameter value to the dictionary
                self.parameters[self.free_parameter_labels[i]].append(value)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the genetic algorithm engine
        self.write_engine()

        # Write the state of the random number generator
        self.write_prng()

    # -----------------------------------------------------------------

    def write_engine(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the state of the genetic engine ...")

        # Save the genetic algorithm
        self.engine.saveto(self.genetic_engine_path_for_generation(self.config.generation_name))

    # -----------------------------------------------------------------

    def write_prng(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the state of the random number generator ...")

        # Save the state of the random generator
        save_state(self.prng_path_for_generation(self.config.generation_name))

# -----------------------------------------------------------------
