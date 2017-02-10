#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.evolve.optimize.stepwise Contains the StepWiseOptimizer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ..core.engine import GeneticEngine
from ...core.tools import filesystem as fs
from ...core.tools.random import save_state, load_state
from ...core.basics.configuration import Configuration
from ..core.adapters import DBFileCSV, DBSQLite
from .optimizer import Optimizer

# -----------------------------------------------------------------

class StepWiseOptimizer(Optimizer):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        This function ...
        :param config:
        """

        # Call the constructor of the base class
        super(StepWiseOptimizer, self).__init__(config)

        # The current population
        self.population = None

        # The scores of the past generation
        self.scores = None

    # -----------------------------------------------------------------

    @classmethod
    def from_directory(cls, dir_path):

        """
        This function ...
        :param dir_path:
        :return:
        """

        # Determine the path to the engine
        engine_path = fs.join(dir_path, "engine.pickle")

        # Load the engine
        engine = GeneticEngine.from_file(engine_path)

        # Determine the path to the random generator
        prng_path = fs.join(dir_path, "prng.pickle")

        # Load the random state
        load_state(prng_path)

        # Determine the path to the configuration
        config_path = fs.join(dir_path, "optimizer.cfg")

        # Load the configuration
        config = Configuration.from_file(config_path)

        # Create the optimizer instance
        optimizer = cls(config)

        # Set the engine
        optimizer.engine = engine

        # Set the database (opening is done during initialization or evolution)
        if fs.extension_of(dir_path, "database") == "csv":
            database_path = fs.join(dir_path, "database.csv")
            optimizer.database = DBFileCSV(filename=database_path, reset=False)
        elif fs.extension_of(dir_path, "database") == "db":
            database_path = fs.join(dir_path, "database.db")
            optimizer.database = DBSQLite(dbname=database_path)
        else: raise ValueError("Unknown extension for database file")

        # Set the path
        #optimizer.path = dir_path
        optimizer.config.output = dir_path

        # Return the optimizer
        return optimizer

    # -----------------------------------------------------------------

    @property
    def initialized(self):

        """
        This function ...
        :return:
        """

        return self.engine is not None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Initialize
        if not self.initialized: self.initialize(**kwargs)

        # 3. Finish
        elif self.config.finish: self.finish()

        # 3. Evolve
        else: self.evolve()

        # 4. Show
        if self.config.show: self.show()

        # 5. Write
        if self.config.write: self.write()

        # 6. Plot
        if self.config.plot: self.plot()

        # Return the current population
        return self.population

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(StepWiseOptimizer, self).setup(**kwargs)

        # Get the engine
        if "engine" in kwargs: self.engine = kwargs.pop("engine")

        # Get the scores
        if "scores" in kwargs: self.scores = kwargs.pop("scores")

        # Get the output path
        if "output" in kwargs: self.config.output = kwargs.pop("output")

        # Set best to None: only when finish_evoluation is run, best should be set
        self.best = None

    # -----------------------------------------------------------------

    def initialize(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Initializing ...")

        # 1. Initialize databse
        self.initialize_database()

        # 2. Initialize genome
        self.initialize_genome(**kwargs)

        # 3. Initialize engine
        self.initialize_engine(**kwargs)

        # 4. Initialize the evolution (the initial generation)
        self.initialize_evolution()

    # -----------------------------------------------------------------

    def initialize_evolution(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Initialize the evolution ...")

        # Initialize the genetic algorithm
        self.engine.initialize()

        # Get the initial population
        self.population = self.engine.get_population()

    # -----------------------------------------------------------------

    def finish(self):

        """
        This function ...
        :return:
        """

        # Set the database adapter again
        self.database.open(self.engine)
        self.engine.setDBAdapter(self.database)

        # Set the scores from the previous generation
        self.set_scores()

        # Set the population
        self.population = self.engine.get_population()

        # Get the best individual
        self.best = self.engine.finish_evolution()

    # -----------------------------------------------------------------

    def evolve(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Evolve ...")

        # Set the database adapter again
        self.database.open(self.engine)
        self.engine.setDBAdapter(self.database)

        # Make sure we don't stop unexpectedly, always increment the number of generations we want
        # (it doesn't matter what the value is)
        self.engine.setGenerations(self.engine.getGenerations() + 1)

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

        # Check if the scores are set
        if self.scores is None: raise ValueError("The scores are not set")

        # Set the scores
        self.engine.set_scores(self.scores)

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

        # Get the new population
        self.population = self.engine.new_population

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

        # Save the engine
        self.write_engine()

        # Save the state of the prng
        self.write_random_state()

        # Save the configuration
        self.write_config()

        # Save the database
        self.write_database()

        # Write the new generation
        self.write_population()

        # Write the best individual
        self.write_best()

    # -----------------------------------------------------------------

    def write_engine(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the engine ...")

        # Determine the path to the engine
        engine_path = self.output_path_file("engine.pickle")

        # Save the engine
        self.engine.saveto(engine_path)

    # -----------------------------------------------------------------

    def write_random_state(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the state of the random number generator ...")

        # Determine the path to the prng
        prng_path = self.output_path_file("prng.pickle")

        # Save the prng state
        save_state(prng_path)

    # -----------------------------------------------------------------

    def write_config(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the configuration file ...")

        # Determine the path
        path = self.output_path_file("optimizer.cfg")

        # Save the configuration
        self.config.saveto(path)

    # -----------------------------------------------------------------

    def write_database(self):

        """
        This function ...
        :param dir_path:
        :return:
        """

        # Inform the user
        log.info("Writing the database ...")

        # Commit all changes and close the database
        self.database.commit_and_close()

    # -----------------------------------------------------------------

    def write_population(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the population genomes ...")

        # Loop over the individuals
        for individual in self.population:

            pass

    # -----------------------------------------------------------------

    def write_best(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the best individual ...")

        pass

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

# -----------------------------------------------------------------
