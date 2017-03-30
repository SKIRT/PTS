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

    def __init__(self, config=None, interactive=False):

        """
        This function ...
        :param config:
        :param interactive:
        """

        # Call the constructor of the base class
        super(StepWiseOptimizer, self).__init__(config, interactive)

        # The current population
        self.population = None

        # The scores of the past generation
        self.scores = None

        # A check for the parameters of the models that are scored
        self.scores_check = None

    # -----------------------------------------------------------------

    @classmethod
    def from_directory(cls, dir_path, run_id):

        """
        This function ...
        :param dir_path:
        :param run_id:
        :return:
        """

        # Debugging
        log.debug("Loading optimizer for run '" + run_id + "' from '" + dir_path + "' ...")

        # Determine the path to the engine
        engine_path = fs.join(dir_path, "engine.pickle")

        # Determine the path to the random generator
        prng_path = fs.join(dir_path, "prng.pickle")

        # Determine the path to the configuration
        config_path = fs.join(dir_path, "optimizer.cfg")

        # Load the statistics (opening is done during initialization or evolution)
        statistics_path = fs.join(dir_path, "statistics.csv")

        # Load the database (opening is done during initialization or evolution)
        database_path = fs.join(dir_path, "database.db")

        # Create the optimizer instance and return it
        return cls.from_paths(dir_path, engine_path, prng_path, config_path, statistics_path, database_path, run_id)

    # -----------------------------------------------------------------

    @classmethod
    def from_paths(cls, output_path, engine_path, prng_path, config_path, statistics_path, database_path, run_id):

        """
        This function ...
        :param output_path:
        :param engine_path:
        :param prng_path:
        :param config_path:
        :param statistics_path:
        :param database_path:
        :param run_id:
        :return:
        """

        # Debugging
        log.debug("")
        log.debug("Loading optimizer for run '" + run_id + "' from:")
        log.debug("")
        log.debug(" - Engine: " + engine_path)
        log.debug(" - PRNG: " + prng_path)
        log.debug(" - Configuration: " + config_path)
        log.debug(" - Statistics: " + statistics_path)
        log.debug(" - Database: " + database_path)
        log.debug(" - Output: " + output_path)
        log.debug("")

        # Load the engine
        engine = GeneticEngine.from_file(engine_path)

        # Load the random state
        load_state(prng_path)

        # Load the configuration
        config = Configuration.from_file(config_path)

        # Create the optimizer instance
        optimizer = cls(config)

        # Set the engine
        optimizer.engine = engine

        # Load the statistics (opening is done during initialization or evolution)
        if not fs.is_file(statistics_path): raise IOError("The statistics file could not be found at '" + statistics_path + "'")
        log.debug("Loading the statistics from '" + statistics_path + "' ...")
        optimizer.statistics = DBFileCSV(filename=statistics_path, reset=False, identify=run_id)

        # Load the database (opening is done during initialization or evolution)
        if not fs.is_file(database_path): raise IOError("The database could not be found at '" + database_path + "'")
        log.debug("Loading the database from '" + database_path + "' ...")
        optimizer.database = DBSQLite(dbname=database_path, resetDB=False, identify=run_id, resetIdentify=False)

        # Set the path
        optimizer.config.output = output_path

        # Set the individual paths
        optimizer.config.writing.engine_path = engine_path
        optimizer.config.writing.prng_path = prng_path
        optimizer.config.writing.config_path = config_path
        optimizer.config.writing.statistics_path = statistics_path
        optimizer.config.writing.database_path = database_path

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

        # Get the scores check
        if "scores_check" in kwargs: self.scores_check = kwargs.pop("scores_check")

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

        # 1. Initialize the statistics table
        self.initialize_statistics()

        # 2. Initialize databse
        self.initialize_database()

        # 3. Initialize genome
        self.initialize_genome(**kwargs)

        # 4. Initialize engine
        self.initialize_engine(**kwargs)

        # 5. Initialize the evolution (the initial generation)
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

        # Set the database and statistics adapters again
        self.set_engine_database_and_statistics()

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
        self.engine.set_scores(self.scores, self.scores_check)

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
        if self.config.writing.engine_path is not None: engine_path = fs.absolute_or_in(self.config.writing.engine_path, self.output_path)
        else: engine_path = self.output_path_file("engine.pickle")

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
        if self.config.writing.prng_path is not None: prng_path = fs.absolute_or_in(self.config.writing.prng_path, self.output_path)
        else: prng_path = self.output_path_file("prng.pickle")

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
        if self.config.writing.config_path is not None: path = fs.absolute_or_in(self.config.writing.config_path, self.output_path)
        else: path = self.output_path_file("optimizer.cfg")

        # Save the configuration
        self.config.saveto(path)

    # -----------------------------------------------------------------

    def write_database(self):

        """
        This function ...
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
