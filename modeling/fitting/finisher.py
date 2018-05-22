#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.finisher Contains the ExplorationFinisher class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.basics.log import log
from ...evolve.optimize.stepwise import StepWiseOptimizer
from .modelgenerators.genetic import set_optimizer_settings, get_last_generation_scores_names_and_check, get_last_generation_name
from ..fitting.sedfitting import SEDFitter
from ..fitting.run import get_ngenerations, has_unfinished_generations, has_unevaluated_generations
from ...core.tools import filesystem as fs
from ...evolve.optimize.stepwise import load_population
from .modelgenerators.genetic import statistics_name, database_name, populations_name, frequency, commit_frequency
from ...evolve.optimize.tables import RecurrenceTable
from pts.core.tools.utils import lazyproperty

# -----------------------------------------------------------------

class ExplorationFinisher(FittingComponent):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ExplorationFinisher, self).__init__(*args, **kwargs)

        # The fitting run
        self.fitting_run = None

        # The scores (only if this is not the initial generation)
        self.scores = None
        self.scores_names = None
        self.scores_check = None

        # The SED fitter
        self.fitter = None

        # The optimizer
        self.optimizer = None

        # The scales for the different parameters
        self.scales = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs: 
        :return: 
        """

        # 2. Evaluate the last generation
        self.evaluate()

        # 3. Do the SED fitting step
        self.fit_sed()

        # Get the scores
        self.get_scores()

        # 4. Set the scores to the optimizer
        self.finish_optimizer()

        # 5. Show
        self.show()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs: 
        :return: 
        """

        # Call the setup function of the base class
        super(ExplorationFinisher, self).setup(**kwargs)

        # Load the fitting run
        self.fitting_run = self.load_fitting_run(self.config.name)

        # Get the scales
        if "scales" in kwargs: self.scales = kwargs.pop("scales")

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

    def evaluate(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Evaluating last generation ...")

        # Check the current number of generations
        current_ngenerations = get_ngenerations(self.config.path, self.config.name)
        # if current_ngenerations <= 1: raise RuntimeError("Need at least one generation after the initial generation to finish the fitting")
        if current_ngenerations == 0: raise RuntimeError("There are no generations")

        # Check if there are unfinished generations
        has_unfinished = has_unfinished_generations(self.config.path, self.config.name)
        if has_unfinished: log.warning("There are unfinished generations, but evaluating finished simulations anyway ...")

        # Check if there are unevaluated generations
        if not has_unevaluated_generations(self.config.path, self.config.name): log.success("All generations have already been evaluated")

    # -----------------------------------------------------------------

    def fit_sed(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Fitting the SED to the finished generations ...")

        # Configuration settings
        config = dict()
        config["name"] = self.config.name

        # Create the SED fitter
        self.fitter = SEDFitter(config)

        # Run the fitter
        self.fitter.run()

        # Success
        if not self.fitting_run.generations_table.has_unfinished: log.success("Succesfully evaluated all generations")

    # -----------------------------------------------------------------

    def get_scores(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Getting the scores ...")

        # Get the scores (set or_initial to False, because this should not happen)
        self.scores, self.scores_names, self.scores_check = get_last_generation_scores_names_and_check(self.fitting_run, or_initial=False)

    # -----------------------------------------------------------------

    def finish_optimizer(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Finishing the optimizer ...")

        # Load the optimizer
        self.load_optimizer()

        # Set settings
        self.set_optimizer_settings()

        # Set path to elitism table for the last generation
        self.optimizer.config.writing.elitism_table_path = self.fitting_run.last_genetic_or_initial_generation.elitism_table_path

        # Set path to scores table for the last generation
        self.optimizer.config.writing.scores_table_path = self.fitting_run.last_genetic_or_initial_generation.scores_table_path

        # Set finish flag
        self.optimizer.config.finish = True

        # Run the optimizer
        self.run_optimizer()

    # -----------------------------------------------------------------

    def load_optimizer(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Loading the optimizer ...")

        # Load the optimizer from files
        self.optimizer = StepWiseOptimizer.from_paths(self.fitting_run.path, self.fitting_run.main_engine_path,
                                                      self.fitting_run.main_prng_path,
                                                      self.fitting_run.optimizer_config_path,
                                                      self.statistics_path,
                                                      self.database_path, self.populations_path, self.config.name,
                                                      statistics_name=statistics_name, database_name=database_name,
                                                      populations_name=populations_name, frequency=frequency,
                                                      commit_frequency=commit_frequency)

    # -----------------------------------------------------------------

    def set_optimizer_settings(self):

        """
        This fucntion ...
        :return: 
        """

        # Inform the user
        log.info("Setting the optimizer settings ...")

        # Set the optimizer settings
        set_optimizer_settings(self.optimizer, self.fitting_run)

    # -----------------------------------------------------------------

    def run_optimizer(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Generating the new models ...")

        evaluator = None
        evaluator_kwargs = None

        # Get the last generation name
        generation_name = get_last_generation_name(self.fitting_run, or_initial=False)

        # Get the parameter minima and maxima
        parameter_minima = self.fitting_run.parameter_minima_for_generation_scalar(generation_name)
        parameter_maxima = self.fitting_run.parameter_maxima_for_generation_scalar(generation_name)

        # Load the previous population
        previous_generation_path = self.fitting_run.last_genetic_or_initial_generation_path
        #previous_population_path = fs.join(previous_generation_path, "population.dat")
        #previous_population = load_population(previous_population_path)
        previous_newborns_path = fs.join(previous_generation_path, "newborns.dat")
        previous_population = load_population(previous_newborns_path)

        # Load the previous recurrence data
        previous_recurrent_path = fs.join(self.fitting_run.last_genetic_generation_path, "recurrence.dat")
        previous_recurrence = RecurrenceTable.from_file(previous_recurrent_path)

        # Run the optimizer
        self.optimizer.run(scores=self.scores, scores_names=self.scores_names, scores_check=self.scores_check, minima=parameter_minima,
                           maxima=parameter_maxima, evaluator=evaluator, evaluator_kwargs=evaluator_kwargs,
                           previous_population=previous_population, previous_recurrence=previous_recurrence,
                           ndigits=self.fitting_run.ndigits_list, nbits=self.fitting_run.nbits_list, scales=self.parameter_scale_list)

    # -----------------------------------------------------------------

    @property
    def population(self):

        """
        This function ...
        :return:
        """

        # The (last) population
        return self.optimizer.population

    # -----------------------------------------------------------------

    @property
    def best(self):

        """
        This function ...
        :return:
        """

        # The best individual
        return self.optimizer.best

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Showing ...")

        # Debugging
        log.debug("Best individual:")
        if log.is_debug: print(self.best)

# -----------------------------------------------------------------
