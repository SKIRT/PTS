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
from ...core.tools.logging import log
from ...evolve.optimize.stepwise import StepWiseOptimizer
from .modelgenerators.genetic import set_optimizer_settings, get_last_generation_scores, get_last_generation_name
from ..fitting.sedfitting import SEDFitter
from ..fitting.run import get_ngenerations, has_unfinished_generations, has_unevaluated_generations

# -----------------------------------------------------------------

class ExplorationFinisher(FittingComponent):

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
        super(ExplorationFinisher, self).__init__(config, interactive)

        # The fitting run
        self.fitting_run = None

        # The scores (only if this is not the initial generation)
        self.scores = None
        self.scores_check = None

        # The SED fitter
        self.fitter = None

        # The optimizer
        self.optimizer = None

        # The (last) population
        self.population = None

        # The best individual
        self.best = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs: 
        :return: 
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Evaluate the last generation
        self.evaluate()

        # 3. Do the SED fitting step
        self.fit_sed()

        # Get the scores
        self.get_scores()

        # 4. Set the scores to the optimizer
        self.finish_optimizer()

        # 5. Get the best individual
        self.get_best()

        # 6. Show
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
        self.scores, self.scores_check = get_last_generation_scores(self.fitting_run, or_initial=False)

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
                                                      self.database_path, self.config.name)

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

        # Run the optimizer
        self.optimizer.run(scores=self.scores, scores_check=self.scores_check, minima=parameter_minima,
                           maxima=parameter_maxima, evaluator=evaluator, evaluator_kwargs=evaluator_kwargs)

    # -----------------------------------------------------------------

    def get_best(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Getting the best individual ...")

        # Set the population
        self.population = self.optimizer.population

        # Get the best individual
        self.best = self.optimizer.best

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
        if log.is_debug(): print(self.best)

# -----------------------------------------------------------------
