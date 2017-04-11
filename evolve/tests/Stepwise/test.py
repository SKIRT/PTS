#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.do.commandline import Command
from pts.evolve.core import reference
from pts.evolve.optimize.optimizer import show_best
from pts.evolve.core.crossovers import G1DListCrossoverOX
from pts.evolve.core.mutators import G1DListMutatorSwap
from pts.core.basics.range import RealRange
from pts.core.test.implementation import TestImplementation
from pts.core.tools.logging import log
from pts.core.tools.loops import repeat

# -----------------------------------------------------------------

description = "finding the maximum of the function defined by Charbonneau (1995) using the StepWiseOptimizer"

# -----------------------------------------------------------------

# Define properties
#nparameters = 20
nparameters = 2
nindividuals = 80
parameter_range = RealRange(0., 1.)
#best_raw_score = float('inf')
best_raw_score = 100
#round_decimal = None
ngenerations = 1000
mutation_rate = 0.03
crossover_rate = 1.0
stats_freq = 100
#mutation_method = "range" # or gaussian, or binary
min_or_max = "maximize"

# -----------------------------------------------------------------

class StepWiseTest(TestImplementation):

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
        super(StepWiseTest, self).__init__(config, interactive)

        # The optimizer
        self.optimizer = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Optimize
        self.optimize()

        # Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(StepWiseTest, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def optimize(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Optimizing ...")

        # Start: launch the initial generation
        self.start()

        # Advance: launch generations 0 -> (n-1)
        repeat(self.advance, self.config.ngenerations)

        # Finish
        self.finish()

    # -----------------------------------------------------------------

    def start(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Starting with a randomly created generation ...")

        # Explore
        self.explore()

    # -----------------------------------------------------------------

    def advance(self):

        """
        THis function ...
        :return: 
        """

        # Inform the user
        log.info("Advancing the fitting with a new generation ...")

        # Load the generations table
        generations = get_generations_table(self.modeling_path, self.fitting_run_name)

        # Check whether there is a generation preceeding this one
        if generations.last_generation_name is None: raise RuntimeError("Preceeding generation cannot be found")

        # Debugging
        log.debug("Previous generation: " + generations.last_generation_name)

        # If some generations have not finished, check the status of and retrieve simulations
        if generations.has_unfinished and self.has_configured_fitting_host_ids: self.synchronize()

        # Debugging
        if generations.has_finished: log.debug(
            "There are finished generations: " + stringify.stringify(generations.finished_generations)[1])
        if has_unevaluated_generations(self.modeling_path, self.fitting_run_name): log.debug(
            "There are unevaluated generations: " +
            stringify.stringify(get_unevaluated_generations(self.modeling_path, self.fitting_run_name))[1])

        # If some generations have finished, fit the SED
        if generations.has_finished and has_unevaluated_generations(self.modeling_path, self.fitting_run_name): self.score() #self.fit_sed()

        # If all generations have finished, explore new generation of models
        if generations.all_finished: self.explore()

    # -----------------------------------------------------------------

    def finish(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Evaluating last generation ...")

        # Check the current number of generations
        current_ngenerations = get_ngenerations(self.modeling_path, self.fitting_run_name)
        # if current_ngenerations <= 1: raise RuntimeError("Need at least one generation after the initial generation to finish the fitting")
        if current_ngenerations == 0: raise RuntimeError("There are no generations")

        # Check if there are unfinished generations
        has_unfinished = has_unfinished_generations(self.modeling_path, self.fitting_run_name)
        if has_unfinished: log.warning(
            "There are unfinished generations, but evaluating finished simulations anyway ...")

        # Check if there are unevaluated generations
        if not has_unevaluated_generations(self.modeling_path, self.fitting_run_name): log.success(
            "All generations have already been evaluated")

        # Do the SED fitting step
        #self.fit_sed()
        self.score()

        # Success
        if not has_unfinished: log.success("Succesfully evaluated all generations")

    # -----------------------------------------------------------------

    def explore(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Exploring the parameter space ...")

        # Configuration settings
        config = dict()
        config["name"] = self.fitting_run_name

        # Create the parameter explorer
        #self.explorer = ParameterExplorer(config)

        # Add an entry to the history
        #self.history.add_entry(ParameterExplorer.command_name())

        # Set the working directory
        #self.explorer.config.path = self.modeling_path

        # Set the remote host IDs
        #self.explorer.config.remotes = self.moderator.host_ids_for_ensemble("fitting")
        #self.explorer.config.attached = self.config.fitting_attached

        # Set the number of generations
        # if self.config.ngenerations is not None: explorer.config.ngenerations = self.config.ngenerations
        # NO: THIS ALWAYS HAVE TO BE ONE: BECAUSE HERE IN THIS CLASS WE ALREADY USE REPEAT(SELF.ADVANCE)
        # IF NGENERATIONS > 1, THE CONTINUOUSOPTIMIZER IS USED INSTEAD OF THE STEPWISEOPTIMIZER
        self.explorer.config.ngenerations = 1

        # Set the number of simulations per generation
        if self.config.nsimulations is not None: self.explorer.config.nsimulations = self.config.nsimulations

        # Set other settings
        #self.explorer.config.npackages_factor = self.config.npackages_factor
        #self.explorer.config.increase_npackages = self.config.increase_npackages
        # explorer.config.refine_wavelengths = self.config.refine_wavelengths
        #self.explorer.config.refine_spectral = self.config.refine_spectral
        # explorer.config.refine_dust = self.config.refine_dust
        #self.explorer.config.refine_spatial = self.config.refine_spatial
        #self.explorer.config.selfabsorption = self.config.selfabsorption
        #self.explorer.config.transient_heating = self.config.transient_heating

        # Set the input
        input_dict = dict()
        if self.parameter_ranges is not None: input_dict["ranges"] = self.parameter_ranges

        # Run the parameter explorer
        self.explorer.run(**input_dict)

        # Mark the end and save the history file
        #self.history.mark_end()
        #self.history.save()

        #def generate_models(self):

        # Inform the user
        log.info("Generating the model parameters ...")

        # Set generator options
        #self.generator.config.ngenerations = self.config.ngenerations
        #self.generator.config.nmodels = self.config.nsimulations

        # Run the model generator
        #self.generator.run(fitting_run=self.fitting_run)

        self.generate()

        self.evaluate()

        self.write_exploration()

    # -----------------------------------------------------------------

    def generate(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Generating the parameter values ...")

        # Reset
        self.optimizer = None

        # Setup
        self.generate_setup()

        # Set parameters
        self.set_parameters()

        # Set settings
        self.set_optimizer_settings()

        # Run the optimizer
        self.run_optimizer()

    # -----------------------------------------------------------------

    def generate_setup(self):

        """
        This function ...
        :return: 
        """

        # Re-invoke existing optimizer run
        if fs.is_file(self.fitting_run.main_engine_path):

            # Load the optimizer from files
            self.optimizer = StepWiseOptimizer.from_paths(self.fitting_run.path,
                                                          self.fitting_run.main_engine_path,
                                                          self.fitting_run.main_prng_path,
                                                          self.fitting_run.optimizer_config_path,
                                                          self.statistics_path, self.database_path,
                                                          self.fitting_run.name)
            # Set initial flag
            self.initial = False

        # New optimizer run
        else:

            # Create a new optimizer and set paths
            self.optimizer = StepWiseOptimizer()
            self.optimizer.config.output = self.fitting_run.path
            self.optimizer.config.writing.engine_path = self.fitting_run.main_engine_path
            self.optimizer.config.writing.prng_path = self.fitting_run.main_prng_path
            self.optimizer.config.writing.config_path = self.fitting_run.optimizer_config_path
            self.optimizer.config.writing.statistics_path = self.statistics_path
            self.optimizer.config.writing.database_path = self.database_path
            self.optimizer.config.run_id = self.fitting_run.name

            # Set initial flag
            self.initial = True

    # -----------------------------------------------------------------

    def set_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the parameter ranges ...")

        # Loop over the free parameters
        for label in self.fitting_run.free_parameter_labels:

            # Use: self.parameter_ranges_for_generation(self, generation_name) ????

            # Get range
            parameter_range = self.fitting_run.free_parameter_ranges[label]

            # Add the parameter
            self.add_parameter(label, parameter_range)

    # -----------------------------------------------------------------

    def set_optimizer_settings(self):

        """
        This function ...
        :return:
        """

        ## In order of optimizer configuration

        # Parameters
        self.optimizer.config.nparameters = self.fitting_run.nfree_parameters

        # Number of generations and the number of individuals per generation
        self.optimizer.config.ngenerations = self.config.ngenerations
        self.optimizer.config.nindividuals = self.config.nmodels

        # User
        self.optimizer.config.mutation_rate = self.fitting_run.genetic_settings.mutation_rate
        self.optimizer.config.crossover_rate = self.fitting_run.genetic_settings.crossover_rate

        # Fixed
        self.optimizer.config.stats_freq = 1
        self.optimizer.config.best_raw_score = 0.

        # User
        self.optimizer.config.round_decimal = self.fitting_run.genetic_settings.round_decimal
        self.optimizer.config.mutation_method = self.fitting_run.genetic_settings.mutation_method

        # User, scaling
        self.optimizer.config.scaling_method = self.fitting_run.genetic_settings.scaling_method

        # User, selector
        self.optimizer.config.selector_method = self.fitting_run.genetic_settings.selector_method

        # Fixed
        self.optimizer.config.min_or_max = "minimize"
        #self.optimizer.config.run_id = self.fitting_run.name # THIS IS NOW DONE IN THE SETUP
        self.optimizer.config.database_frequency = 1
        self.optimizer.config.statistics_frequency = 1

        # Fixed
        #self.optimizer.config.output = self.fitting_run.path

        # Fixed
        self.optimizer.config.elitism = True

        # Fixed
        self.optimizer.config.nelite_individuals = self.fitting_run.genetic_settings.nelite_individuals

        # Set heterogeneous flag
        self.optimizer.config.heterogeneous = True

    # -----------------------------------------------------------------

    def run_optimizer(self):

        # Inform the user
        log.info("Generating the new models ...")

        # Set the scores
        if not self.initial: self.set_scores()

        # Run the optimizer
        self.optimizer.run(scores=self.scores, scores_check=self.scores_check, minima=self.parameter_minima_scalar,
                           maxima=self.parameter_maxima_scalar, evaluator=self.evaluator,
                           evaluator_kwargs=self.evaluator_kwargs)

        # Get the parameter values of the new models
        self.get_model_parameters()

    # -----------------------------------------------------------------

    def set_scores(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting scores from previous generation ...")

        # Load the parameters table from the previous generation
        parameters_table = self.fitting_run.parameters_table_for_generation(self.fitting_run.last_genetic_or_initial_generation_name)

        # Load the chi squared table from the previous generation
        chi_squared_table = self.fitting_run.chi_squared_table_for_generation(self.fitting_run.last_genetic_or_initial_generation_name)

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
        for label in self.fitting_run.free_parameter_labels:
            values = parameters_table[label]
            check.append(values)

        # Set the scores
        self.scores = scores
        self.scores_check = check

    # -----------------------------------------------------------------

    def get_model_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the model parameters ...")

        # Loop over the individuals of the population
        for individual in self.optimizer.population:

            # Loop over all the genes (parameters)
            for i in range(len(individual)):

                # Get the parameter value
                value = individual[i]

                # Add the parameter value to the dictionary
                self.parameters[self.fitting_run.free_parameter_labels[i]].append(value)

    # -----------------------------------------------------------------

    def evaluate(self):

        """
        THis function ...
        :return: 
        """

        # Inform the user
        log.info("Evaluating ...")

        nmodels = self.generator.nmodels

        # Loop over the different parameter combinations
        for i in range(self.nmodels):

            # Get the parameter values as a dictionary
            parameter_values = get_parameter_values_from_generator(self.generator.parameters, i, self.fitting_run)

    # -----------------------------------------------------------------

    def write_exploration(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Writing after exploration step ...")

        self.write_generation()

        self.write_parameters()

        self.write_chi_squared()

    # -----------------------------------------------------------------

    def write_generation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing generation info ...")

        # Add an entry to the generations table
        self.fitting_run.generations_table.add_entry(self.generation, self.ranges)

        # Save the table
        self.fitting_run.generations_table.save()

    # -----------------------------------------------------------------

    def write_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the model parameters table ...")

        # Save the parameters table
        self.parameters_table.saveto(self.generation.parameters_table_path)

    # -----------------------------------------------------------------

    def write_chi_squared(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the chi squared table ...")

        # Save the chi squared table
        self.chi_squared_table.saveto(self.generation.chi_squared_table_path)

    # -----------------------------------------------------------------

    def score(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Scoring engine ...")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Writing ...")

# -----------------------------------------------------------------

def eval_func_xy(x, y):

    """
    This function ...
    :param x:
    :param y:
    :return:
    """

    # Calculate z and return
    z = (16 * x * (1. - x) * y * (1. - y) * np.sin(2. * np.pi * x) * np.sin(2. * np.pi * y) )**2
    return z

# -----------------------------------------------------------------

def eval_func(chromosome, **kwargs):

    """
    The evaluation function
    """

    # Get x and y
    x = chromosome[0]
    y = chromosome[1]

    # Return the value of the function
    return eval_func_xy(x, y)

# -----------------------------------------------------------------

def evolve_callback(ga_engine, **kwargs):

    """
    This function ...
    :param ga_engine:
    :param kwargs:
    :return:
    """

    # Get the coords
    #coords = kwargs.pop("coordinates")

    #if ga_engine.currentGeneration % 10 == 0:
    #    best = ga_engine.bestIndividual()
    #    write_tour_to_img(coords, best, "tsp_result_%d.png" % (ga_engine.currentGeneration,))
    #return False


    #return coords, cm

    pass

# -----------------------------------------------------------------
# SETUP FUNCTION
# -----------------------------------------------------------------

def setup(temp_path):

    """
    This function ...
    :param temp_path:
    """

    return

# -----------------------------------------------------------------
# OPTIMIZE
# -----------------------------------------------------------------

# Settings
settings_optimize = dict()
settings_optimize["output"] = None
settings_optimize["nparameters"] = nparameters
settings_optimize["nindividuals"] = nindividuals
settings_optimize["parameter_range"] = parameter_range
settings_optimize["best_raw_score"] = best_raw_score
#settings_optimize["round_decimal"] = round_decimal
settings_optimize["ngenerations"] = ngenerations
settings_optimize["mutation_rate"] = mutation_rate
settings_optimize["crossover_rate"] = crossover_rate
settings_optimize["stats_freq"] = stats_freq
#settings_optimize["mutation_method"] = mutation_method
settings_optimize["min_or_max"] = min_or_max

# Input
input_optimize = dict()
#input_optimize["genome"] = genome
input_optimize["evaluator"] = eval_func
#input_optimize["initializator"] = G1DListTSPInitializator
input_optimize["mutator"] = G1DListMutatorSwap
input_optimize["crossover"] = G1DListCrossoverOX
input_optimize["callback"] = evolve_callback
#input_optimize["adapter"] = sqlite_adapter

# Create dictionary for extra arguments to the evalutor function
#input_optimize["evaluator_kwargs"] = {"distances": cm}
#input_optimize["callback_kwargs"] = {"coordinates": coords}

# -----------------------------------------------------------------
# TEST FUNCTION
# -----------------------------------------------------------------

def test(temp_path):

    """
    This function ...
    """

    # Remove 'callback' from the settings dictionary: we don't want to plot now
    #if "callback" in input_optimize: del input_optimize["callback"]

    # Solve the problem with the original Pyevolve implementation
    best = reference.call(settings_optimize, input_optimize)

    # Show the best individual
    show_best(best)

# -----------------------------------------------------------------
