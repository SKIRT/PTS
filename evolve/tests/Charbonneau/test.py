#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import inspect
import numpy as np

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.do.commandline import Command
from pts.evolve.core import reference
from pts.evolve.optimize.optimizer import show_best
from pts.evolve.core.crossovers import G1DListCrossoverOX
from pts.evolve.core.mutators import G1DListMutatorSwap
from pts.core.basics.animation import Animation
from pts.core.basics.range import RealRange
from pts.core.test.implementation import TestImplementation
from pts.core.basics.log import log

# -----------------------------------------------------------------

this_path = fs.absolute_path(inspect.stack()[0][1])
this_dir_path = fs.directory_of(this_path)

# -----------------------------------------------------------------

description = "finding the maximum of the function defined by Charbonneau (1995)"

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

class CharbonneauTest(TestImplementation):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(CharbonneauTest, self).__init__(*args, **kwargs)

        # The optimizer
        self.optimizer = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Optimize
        self.optimize()

        # Make animation
        self.make_animation()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(CharbonneauTest, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def optimize(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Optimizing ...")

        # Construct the command
        command = Command("optimize_continuous", "finding the maximum of the function defined by Charbonneau (1995)",
                           settings_optimize, input_optimize, cwd=".")

        # Add the command
        #commands.append(optimize)

        self.optimizer = self.run_command(command)

    # -----------------------------------------------------------------

    def make_animation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making animation ...")

        # Get the best
        best = self.optimizer.best

        # Determine path
        temp_path = self.optimizer.config.path
        filepath = fs.join(temp_path, "best.png")

        # Make plot of best

        # Create animated gif
        animation = Animation()

        # Open the images present in the directory
        for path in fs.files_in_path(temp_path, extension="png", exact_not_name="best"):
            # Add frame to the animation
            animation.add_frame_from_file(path)

        # Save the animation
        animation_path = fs.join(temp_path, "animation.gif")
        animation.saveto(animation_path)

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
