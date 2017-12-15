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
from ....core.basics.log import log
from .generator import ModelGenerator
from ....core.tools import filesystem as fs
from ....evolve.optimize.stepwise import StepWiseOptimizer, load_population
from ....evolve.optimize.continuous import ContinuousOptimizer
from ..evaluate import evaluate
from ....core.tools.stringify import tostr
from ....core.tools import sequences
from ....core.tools import numbers
from ....evolve.optimize.tables import RecurrenceTable
from pts.core.tools.utils import lazyproperty

# -----------------------------------------------------------------

database_name = "database"
statistics_name = "statistics file"
populations_name = "populations file"

frequency = 1
commit_frequency = 1

# -----------------------------------------------------------------

class GeneticModelGenerator(ModelGenerator):

    """
    This function ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param interactive:
        """

        # Call the constructor of the base class
        super(GeneticModelGenerator, self).__init__(*args, **kwargs)

        # The scores (only if this is not the initial generation)
        self.scores = None
        self.scores_names = None
        self.scores_check = None

        # The optimizer
        self.optimizer = None

        # Whether or not this is the initial generation
        self.initial = None

        # The evaluator function and its kwargs
        self.evaluator = None
        self.evaluator_kwargs = None

        # The parameter ranges
        self.parameter_ranges = None

        # A dictionary of fixed parameter values for the initial generation
        self.fixed_initial_parameters = None

        # The parameter values of the intial populaiton
        self.initial_parameters = None

        # Previous population
        self.previous_population = None

        # Previous recurrent
        self.previous_recurrence = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(GeneticModelGenerator, self).setup(**kwargs)

        # If the number of generations in one run is more than one
        if self.config.ngenerations > 1: self.create_continuous_optimizer()

        # Only one generation at a time
        else: self.set_stepwise_optimizer()

        # Get the parameter ranges
        if "parameter_ranges" in kwargs: self.parameter_ranges = kwargs.pop("parameter_ranges")

        # Get the fixed initial parameter values
        if "fixed_initial_parameters" in kwargs: self.fixed_initial_parameters= kwargs.pop("fixed_initial_parameters")

        # Set parameters
        self.set_parameters()

        # Set settings
        self.set_optimizer_settings()

    # -----------------------------------------------------------------

    def create_continuous_optimizer(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Creating an optimizer that runs multiple generations at once ...")

        # Create a new continuous optimizer
        self.optimizer = ContinuousOptimizer()
        self.optimizer.config.output = self.fitting_run.path
        # self.optimizer.config.writing.engine_path = self.fitting_run.main_engine_path
        # self.optimizer.config.writing.prng_path = self.fitting_run.main_prng_path
        # self.optimizer.config.writing.config_path = self.fitting_run.optimizer_config_path
        # self.optimizer.config.writing.statistics_path = self.statistics_path
        # self.optimizer.config.writing.database_path = self.database_path

        self.optimizer.config.run_id = self.fitting_run.name

        # Set initial flag
        self.initial = True

        # Set the evaluator function
        self.evaluator = evaluate

        # Set the evaluator kwargs
        self.evaluator_kwargs = dict()
        self.evaluator_kwargs["fitting_run"] = self.fitting_run  # self.fitting_run.name

    # -----------------------------------------------------------------

    def set_stepwise_optimizer(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Setting the optimizer that runs one generation at a time ...")

        # Re-invoke existing optimizer run
        if fs.is_file(self.fitting_run.main_engine_path): self.load_optimizer()

        # New optimizer run
        else: self.create_optimizer()

        # Set elitism table path -> directory of previous generation!
        if self.initial: elitism_path = None
        else: elitism_path = self.fitting_run.last_genetic_or_initial_generation.elitism_table_path # fs.join(self.fitting_run.last_genetic_or_initial_generation_path, "elitism.dat")

        # Set scores table path -> directory of previous generation!
        if self.initial: scores_path = None
        else: scores_path = self.fitting_run.last_genetic_or_initial_generation.scores_table_path

        # Set generation specific paths
        self.optimizer.config.writing.input_path = self.generation.optimizer_input_path
        self.optimizer.config.writing.newborns_path = self.generation.newborns_path
        self.optimizer.config.writing.parents_path = self.generation.parents_path
        self.optimizer.config.writing.crossover_table_path = self.generation.crossover_table_path
        self.optimizer.config.writing.recurrent_path = self.generation.recurrence_path
        self.optimizer.config.writing.elitism_table_path = elitism_path
        self.optimizer.config.writing.scores_table_path = scores_path

        # NEW: ADD EXTRA PATHS FOR WRITING ENGINE, PRNG AND CONFIG
        self.optimizer.config.writing.engine_path = [self.optimizer.config.writing.engine_path, self.generation.engine_path]
        self.optimizer.config.writing.prng_path = [self.optimizer.config.writing.prng_path, self.generation.prng_path]
        self.optimizer.config.writing.config_path = [self.optimizer.config.writing.config_path, self.generation.optimizer_config_path]

        # Set recurrence settings
        self.optimizer.config.check_recurrence = self.config.check_recurrence
        self.optimizer.config.recurrence_rtol = self.config.recurrence_rtol
        self.optimizer.config.recurrence_atol = self.config.recurrence_atol

    # -----------------------------------------------------------------

    def load_optimizer(self):

        """
        This function ...
        :param self: 
        :return: 
        """

        # Inform the user
        log.info("Loading the optimizer from the previous generation ....")

        # Load the optimizer from files
        self.optimizer = StepWiseOptimizer.from_paths(self.fitting_run.path,
                                                      self.fitting_run.main_engine_path,
                                                      self.fitting_run.main_prng_path,
                                                      self.fitting_run.optimizer_config_path,
                                                      self.statistics_path, self.database_path,
                                                      self.populations_path, self.fitting_run.name,
                                                      statistics_name=statistics_name, database_name=database_name,
                                                      populations_name=populations_name, frequency=frequency,
                                                      commit_frequency=commit_frequency)
        # Set initial flag
        self.initial = False

        # Load the previous population
        previous_generation_path = self.fitting_run.last_genetic_or_initial_generation_path
        #previous_population_path = fs.join(previous_generation_path, "population.dat")
        previous_parents_path = fs.join(previous_generation_path, "parents.dat")
        previous_newborns_path = fs.join(previous_generation_path, "newborns.dat")
        if fs.is_file(previous_newborns_path): self.previous_population = load_population(previous_newborns_path)
        else: self.previous_population = load_population(previous_parents_path)
        #self.previous_population = load_population(previous_population_path)

        # Load the previous recurrent data, BUT NOT WHEN THE PREVIOUS GENERATION WAS THE INTIAL ONE ?
        if not self.fitting_run.last_is_initial:
            previous_recurrent_path = fs.join(self.fitting_run.last_genetic_generation_path, "recurrence.dat")
            self.previous_recurrence = RecurrenceTable.from_file(previous_recurrent_path)

    # -----------------------------------------------------------------

    def create_optimizer(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Creating a new stepwise optimizer ...")

        # Create a new optimizer and set paths
        self.optimizer = StepWiseOptimizer()
        self.optimizer.config.output = self.fitting_run.path

        # Set basic paths
        self.optimizer.config.writing.engine_path = self.fitting_run.main_engine_path
        self.optimizer.config.writing.prng_path = self.fitting_run.main_prng_path
        self.optimizer.config.writing.config_path = self.fitting_run.optimizer_config_path

        # Set other paths
        self.optimizer.config.writing.statistics_path = self.statistics_path
        self.optimizer.config.writing.database_path = self.database_path
        self.optimizer.config.writing.populations_path = self.populations_path

        # Set fitting run ID
        self.optimizer.config.run_id = self.fitting_run.name

        # Set initial flag
        self.initial = True

        # Set previous population to None
        self.previous_population = None

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

            # Get range (=default range)
            default_parameter_range = self.fitting_run.free_parameter_ranges[label]

            # Get the name of the last generation
            # DOESN'T WORK: SEE BELOW WHY
            #generation_name = get_last_generation_name(self.fitting_run)

            # DOESN'T WORK: the generations table has not been updated at this point: only when write() is called from the ParameterExplorer (parent to this class)
            # parameter_range = self.fitting_run.parameter_ranges_for_generation(generation_name, label)

            # Set the parameter range
            if self.parameter_ranges is not None and label in self.parameter_ranges: parameter_range = self.parameter_ranges[label]
            else: parameter_range = default_parameter_range

            # Debugging
            log.debug("Parameter '" + label + "':")
            log.debug("")
            log.debug(" - default range: " + tostr(default_parameter_range, fancy=True))
            log.debug(" - range for this generation: " + tostr(parameter_range, fancy=True))
            log.debug("")

            # Add the parameter
            self.add_parameter(label, parameter_range)

    # -----------------------------------------------------------------

    def set_optimizer_settings(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the settings for the optimizer ...")

        # Set the settings
        set_optimizer_settings(self.optimizer, self.fitting_run, self.config.ngenerations, self.config.nmodels)

    # -----------------------------------------------------------------

    @property
    def generate_initial_manual(self):

        """
        This poret ekgpoege
        :return: 
        """

        if self.fixed_initial_parameters is not None: return True
        else: return self.config.manual_initial_generation

    # -----------------------------------------------------------------

    def generate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the new models ...")

        # Set the scores
        if not self.initial: self.set_scores()

        # Generate parameters of the initial generation
        if self.initial and self.generate_initial_manual: self.generate_initial_parameters()

        # Run the optimizer
        self.optimizer.run(scores=self.scores, scores_check=self.scores_check, minima=self.parameter_minima_scalar,
                           scores_names=self.scores_names,
                           maxima=self.parameter_maxima_scalar, evaluator=self.evaluator,
                           evaluator_kwargs=self.evaluator_kwargs, initial_parameters=self.initial_parameters,
                           previous_population=self.previous_population, previous_recurrence=self.previous_recurrence,
                           ndigits=self.fitting_run.ndigits_list, nbits=self.fitting_run.nbits_list, scales=self.parameter_scale_list)

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

        # Get the scores (and names and check)
        self.scores, self.scores_names, self.scores_check = get_last_generation_scores_names_and_check(self.fitting_run)

    # -----------------------------------------------------------------

    @property
    def nfixed_initial_parameters(self):

        """
        This property ...
        :return: 
        """

        # No fixed initial parameters
        if self.fixed_initial_parameters is None: return 0
        else: return len(self.fixed_initial_parameters[self.fixed_initial_parameters.keys()[0]])

    # -----------------------------------------------------------------

    def generate_initial_parameters(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Generating the parameters of the initial generation ...")

        # Initialize a list to contain the parameter sets
        self.initial_parameters = []

        # Add fixed initial parameter sets if desired
        if self.fixed_initial_parameters is not None: self.generate_initial_parameters_fixed()

        # Add random initial parameter sets
        if self.config.manual_initial_generation_method == "random": self.generate_initial_parameters_random()

        # Add initial parameter sets based on a grid
        elif self.config.manual_initial_generation_method == "grid": self.generate_initial_parameters_grid()

        # Invalid method
        else: raise ValueError("Invalid manual initial generation method: " + self.config.manual_initial_generation_method)

    # -----------------------------------------------------------------

    def generate_initial_parameters_fixed(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Generating parameter sets based on the fixed initial parameters ...")

        # Loop over the number of fixed initial parameter sets
        for i in range(self.nfixed_initial_parameters):

            # Initialize the parameter set
            initial_parameters_model = []

            # Get the value for each free parameters
            for label in self.fitting_run.free_parameter_labels:

                # Get the quantity
                quantity = self.fixed_initial_parameters[label][i]

                # Get only the scalar value
                if label in self.fitting_run.parameter_units and self.fitting_run.parameter_units[label] is not None:
                    unit = self.fitting_run.parameter_units[label]
                    value = quantity.to(unit).value
                else: value = quantity

                # Add the value to the parameter set
                initial_parameters_model.append(value)

            # Add the parameter set
            self.initial_parameters.append(initial_parameters_model)

    # -----------------------------------------------------------------

    @property
    def single_parameter_scale(self):
        
        """
        This function ..
        :return: 
        """

        return self.scales is None

    # -----------------------------------------------------------------

    @property
    def multiple_parameter_scales(self):

        """
        This function ...
        :return: 
        """

        return self.scales is not None

    # -----------------------------------------------------------------

    def scale_for_parameter(self, label):

        """
        This function ...
        :param label: 
        :return: 
        """

        if self.multiple_parameter_scales: return self.scales[label]
        else: return self.config.default_scale

    # -----------------------------------------------------------------

    @lazyproperty
    def parameter_scale_list(self):

        """
        This function ...
        :return: 
        """

        scales = []
        for label in self.fitting_run.free_parameter_labels: scales.append(self.scale_for_parameter(label))
        return scales

    # -----------------------------------------------------------------

    def generate_initial_parameters_random(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Generating initial parameter sets based on random number generation ...")

        # Loop over the number of required models minus the number of fixed model parameter sets
        for _ in range(self.config.nmodels - self.nfixed_initial_parameters):

            # Initialize list for parameter set
            initial_parameters_model = []

            # Set the list values
            for label_index, label in enumerate(self.fitting_run.free_parameter_labels):

                # Get range limits
                range_min = self.parameter_minima_scalar[label_index]
                range_max = self.parameter_maxima_scalar[label_index]

                # Linear scale
                if self.scale_for_parameter(label) == "linear": random = numbers.random_linear(range_min, range_max)

                # Logarithmic scale
                elif self.scale_for_parameter(label) == "logarithmic": random = numbers.random_logarithmic(range_min, range_max)

                # Invalid scale
                else: raise ValueError("Invalid scale: " + str(self.scale_for_parameter(label)))

                # Add value to the parameter set
                initial_parameters_model.append(random)

            # Add the parameter set
            self.initial_parameters.append(initial_parameters_model)

    # -----------------------------------------------------------------

    def generate_initial_parameters_grid(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Generating initial parameter sets based on a grid in the n-d parameter space ...")

        # Generate grid points
        if self.single_parameter_scale: grid_points = self.generate_grid_points_one_scale(self.config.default_scale, most_sampled=self.most_sampled_parameters, weights=self.sampling_weights)
        else: grid_points = self.generate_grid_points_different_scales(self.scales, most_sampled=self.most_sampled_parameters, weights=self.sampling_weights)

        # Convert into lists
        grid_points = self.grid_points_to_lists(grid_points)

        # Create iterator of combinations
        iterator = sequences.iterate_lists_combinations(*grid_points)

        # Generate the initial parameter sets
        # Loop over the number of required models minus the number of fixed model parameter sets
        for index in range(self.config.nmodels - self.nfixed_initial_parameters):

            # The next combination
            initial_parameters_model = list(iterator.next()) # returns tuple

            # Add the parameter set
            self.initial_parameters.append(initial_parameters_model)

    # -----------------------------------------------------------------

    @property
    def individual_names(self):

        """
        This function ...
        :return: 
        """

        return self.optimizer.new_individual_names

    # -----------------------------------------------------------------

    def get_model_parameters(self):

        """
        THis function ...
        :return: 
        """

        # Inform the user
        log.info("Getting the model parameters ...")

        # Loop over the individual names
        for name in self.individual_names:

            # Get the parameters
            for index, parameter_value in enumerate(self.optimizer.get_parameters(name)):

                # Add the parameter value to the dictionary
                self.parameters[self.fitting_run.free_parameter_labels[index]][name] = parameter_value

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the state of the optimizer
        #self.write_optimizer()

# -----------------------------------------------------------------

def set_optimizer_settings(optimizer, fitting_run, ngenerations=None, nmodels=None):

    """
    This function ...
    :param optimizer: 
    :param fitting_run: 
    :param ngenerations:
    :param nmodels:
    :return: 
    """

    ## In order of optimizer configuration

    # NEW: GENOME TYPE
    optimizer.config.genome_type = fitting_run.genetic_settings.genome_type
    optimizer.config.genome_dimension = 1

    # Parameters
    optimizer.config.nparameters = fitting_run.nfree_parameters

    # Number of generations and the number of individuals per generation
    optimizer.config.ngenerations = ngenerations
    optimizer.config.nindividuals = nmodels

    # User
    optimizer.config.mutation_rate = fitting_run.genetic_settings.mutation_rate
    optimizer.config.crossover_rate = fitting_run.genetic_settings.crossover_rate

    # Fixed
    optimizer.config.stats_freq = 1
    optimizer.config.best_raw_score = 0.

    # User
    optimizer.config.round_decimal = fitting_run.genetic_settings.round_decimal
    optimizer.config.mutation_method = fitting_run.genetic_settings.mutation_method
    optimizer.config.binary_mutation_method = fitting_run.genetic_settings.binary_mutation_method

    # User, scaling
    optimizer.config.scaling_method = fitting_run.genetic_settings.scaling_method

    # User, selector
    optimizer.config.selector_method = fitting_run.genetic_settings.selector_method

    # Fixed
    optimizer.config.min_or_max = "minimize"
    # self.optimizer.config.run_id = self.fitting_run.name # THIS IS NOW DONE IN THE SETUP
    optimizer.config.database_frequency = frequency
    optimizer.config.database_commit_frequency = commit_frequency
    optimizer.config.statistics_frequency = frequency
    optimizer.config.populations_frequency = frequency

    # Set names
    optimizer.config.database_name = database_name
    optimizer.config.statistics_name = statistics_name
    optimizer.config.populations_name = populations_name

    # Fixed
    # self.optimizer.config.output = self.fitting_run.path

    # Fixed
    optimizer.config.elitism = True

    # Fixed
    optimizer.config.nelite_individuals = fitting_run.genetic_settings.nelite_individuals

    # Set heterogeneous flag
    optimizer.config.heterogeneous = True

    # Set named_individuals flag
    optimizer.config.named_individuals = True

    # Recurrency checking
    #optimizer.config.check_recurrence = True
    #optimizer.config.recurrence_rtol = 1e-5
    #optimizer.config.recurrence_atol = 1e-8

    # ADVANCED: GRAY CODING
    optimizer.config.gray_code = fitting_run.genetic_settings.gray_code

    # ADVANCED: THE RANDOM SEED
    optimizer.config.seed = fitting_run.random_seed

# -----------------------------------------------------------------

def get_last_generation_name(fitting_run, or_initial=True):

    """
    This function ...
    :param fitting_run: 
    :param or_initial: 
    :return: 
    """

    # Get the corresponding generation name
    if or_initial: generation_name = fitting_run.last_genetic_or_initial_generation_name
    else: generation_name = fitting_run.last_genetic_generation_name

    # Return the name
    return generation_name

# -----------------------------------------------------------------

def get_last_generation_scores_and_check(fitting_run, or_initial=True):

    """
    This function ...
    :param fitting_run:
    :param or_initial:
    :return: 
    """

    # Get the generation name
    generation_name = get_last_generation_name(fitting_run, or_initial=or_initial)

    # Load the parameters table from the previous generation
    parameters_table = fitting_run.parameters_table_for_generation(generation_name)

    # Load the chi squared table from the previous generation
    chi_squared_table = fitting_run.chi_squared_table_for_generation(generation_name)

    # List of chi squared values in the same order as the parameters table
    chi_squared_values = []

    # The checks
    check = []

    # Check whether the chi-squared and parameter tables match
    for i in range(len(parameters_table)):

        simulation_name = parameters_table["Simulation name"][i]
        chi_squared = chi_squared_table.chi_squared_for(simulation_name)
        chi_squared_values.append(chi_squared)

        # Check individual values with parameter table of the last generation
        check_individual = []
        for label in fitting_run.free_parameter_labels:
            value = parameters_table[label][i]
            check_individual.append(value)
        check.append(check_individual)

    # Set the scores
    scores = chi_squared_values

    # Return
    return scores, check

# -----------------------------------------------------------------

def get_last_generation_scores_and_names(fitting_run, or_initial=True):

    """
    This function ...
    :param fitting_run:
    :param or_initial:
    :return:
    """

    # Get the generation name
    generation_name = get_last_generation_name(fitting_run, or_initial=or_initial)

    # Load the parameters table from the previous generation
    parameters_table = fitting_run.parameters_table_for_generation(generation_name)

    # Load the chi squared table from the previous generation
    chi_squared_table = fitting_run.chi_squared_table_for_generation(generation_name)

    # Load the individuals table from the previous generation
    individuals_table = fitting_run.individuals_table_for_generation(generation_name)

    # List of chi squared values in the same ordere as the parameters table
    chi_squared_values = []

    # The names
    names = []

    # Loop over the simulations
    for i in range(len(parameters_table)):

        simulation_name = parameters_table["Simulation name"][i]
        chi_squared = chi_squared_table.chi_squared_for(simulation_name)
        chi_squared_values.append(chi_squared)

        # Get the individual name
        individual_name = individuals_table.get_individual_name(simulation_name)

        # Add the name
        names.append(individual_name)

    # Set the scores
    scores = chi_squared_values

    # Return
    return scores, names

# -----------------------------------------------------------------

def get_last_generation_scores_names_and_check(fitting_run, or_initial=True):

    """
    THis function ...
    :param fitting_run:
    :param or_initial:
    :return:
    """

    # Get the generation name
    generation_name = get_last_generation_name(fitting_run, or_initial=or_initial)

    # Load the parameters table from the previous generation
    parameters_table = fitting_run.parameters_table_for_generation(generation_name)

    # Load the chi squared table from the previous generation
    chi_squared_table = fitting_run.chi_squared_table_for_generation(generation_name)

    # Load the individuals table from the previous generation
    individuals_table = fitting_run.individuals_table_for_generation(generation_name)

    # List of chi squared values in the same ordere as the parameters table
    chi_squared_values = []

    # The names
    names = []

    # The checks
    check = []

    # Loop over the simulations
    for i in range(len(parameters_table)):

        simulation_name = parameters_table["Simulation name"][i]
        chi_squared = chi_squared_table.chi_squared_for(simulation_name)
        chi_squared_values.append(chi_squared)

        # Get the individual name
        individual_name = individuals_table.get_individual_name(simulation_name)

        # Add the name
        names.append(individual_name)

        # Check individual values with parameter table of the last generation
        check_individual = []
        for label in fitting_run.free_parameter_labels:
            value = parameters_table[label][i]
            check_individual.append(value)
        check.append(check_individual)

    # Set the scores
    scores = chi_squared_values

    # Return
    return scores, names, check

# -----------------------------------------------------------------
