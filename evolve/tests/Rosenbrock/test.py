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
from matplotlib import pyplot as plt

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.do.commandline import Command
from pts.core.basics.range import IntegerRange
from pts.evolve.core import reference
from pts.evolve.optimize.optimizer import show_best
from pts.core.test.implementation import TestImplementation

# -----------------------------------------------------------------

this_path = fs.absolute_path(inspect.stack()[0][1])
this_dir_path = fs.directory_of(this_path)

# -----------------------------------------------------------------

description = "optimizing the Rosenbrock function, a deceptive function"

# -----------------------------------------------------------------

# Define properties
nparameters = 20
#nparameters = 2
nindividuals = 80
parameter_range = IntegerRange(0, 10)
best_raw_score = 0.0
round_decimal = None
ngenerations = 4000
mutation_rate = 0.2
crossover_rate = None
stats_freq = 200
mutation_method = "range" # or gaussian, or binary
min_or_max = "minimize"

# -----------------------------------------------------------------

class RosenbrockTest(TestImplementation):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(RosenbrockTest, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Optimize
        self.optimize()

        # Reference
        self.reference()

    # -----------------------------------------------------------------

    def optimize(self):

        """
        This function ...
        :return:
        """

        # Settings
        settings_optimize = dict()
        settings_optimize["output"] = None
        settings_optimize["nparameters"] = nparameters
        settings_optimize["nindividuals"] = nindividuals
        settings_optimize["parameter_range"] = parameter_range
        settings_optimize["best_raw_score"] = best_raw_score
        settings_optimize["round_decimal"] = round_decimal
        settings_optimize["ngenerations"] = ngenerations
        settings_optimize["mutation_rate"] = mutation_rate
        settings_optimize["crossover_rate"] = crossover_rate
        settings_optimize["stats_freq"] = stats_freq
        settings_optimize["mutation_method"] = mutation_method
        settings_optimize["min_or_max"] = min_or_max

        # Input
        input_optimize = dict()
        input_optimize["evaluator"] = rosenbrock

        # -----------------------------------------------------------------

        # Create plot
        # ax = create_plot()

        # Create callback to plot the best individuals for each generation
        input_optimize["callback"] = add_best_to_plot
        # input_optimize["callback_kwargs"] = {"ax": ax}

        # -----------------------------------------------------------------

        # Construct the command
        optimize = Command("optimize", "optimize the Rosenbrock function", settings_optimize, input_optimize, cwd=".", finish=show_plot)

        # Add the command
        #commands.append(optimize)

    # -----------------------------------------------------------------

    def reference(self):

        """
        This function ...
        :return:
        """

        # Solve the problem with the original Pyevolve implementation
        best = reference.call(settings_optimize, input_optimize)

        # Show the best individual
        show_best(best)

# -----------------------------------------------------------------

def rosenbrock(xlist):

    """
    This is the Rosenbrock function, a deceptive function
    :param xlist:
    :return:
    """

    sum1 = 0
    for x in xrange(1, len(xlist)):

      sum1 += 100.0 * (xlist[x] - xlist[x-1]**2)**2 + (1 - xlist[x-1])**2

    # Return the raw score
    return sum1

# -----------------------------------------------------------------

def add_best_to_plot(ga_engine, **kwargs):

    """
    This function is called after each step (generation)
    :param ga_engine:
    :return:
    """

    ax = kwargs.pop("ax")

    best = ga_engine.bestIndividual()

    index = ga_engine.currentGeneration

    if index % 50 != 0: return

    generation_label = "Generation " + str(index)

    x = best.genomeList[0]
    y = best.genomeList[1]

    color = "#000000"

    # Plot the point
    ax.plot(x, y, 'o-', color=color, label=generation_label, alpha=0.8, lw=2, markersize=5, mew=1, mec=color, mfc='none')

# -----------------------------------------------------------------

def create_plot():

    """
    This function ...
    :return:
    """

    _, ax = plt.subplots(1, 1)

    # make a contour plot of the rosenbrock function surface.
    X, Y = np.meshgrid(np.linspace(-1.3, 1.3, 31), np.linspace(-0.9, 1.7, 31))
    Z = 100 * (Y - X ** 2) ** 2 + (1 - X) ** 2
    ax.plot([1], [1], 'x', mew=3, markersize=10, color='#111111')
    ax.contourf(X, Y, Z, np.logspace(-1, 3, 31), cmap='gray_r')

    # Set axes limits
    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-0.9, 1.7)

    # Return the axes
    return ax

# -----------------------------------------------------------------

def show_plot(*args, **kwargs):

    """
    This function ...
    :param args:
    :param kwargs:
    :return:
    """

    # Show
    plt.legend(loc='lower right')
    plt.show()

# -----------------------------------------------------------------
# TEST FUNCTION
# -----------------------------------------------------------------

def test(temp_path):

    """
    This function ...
    """

    pass

# -----------------------------------------------------------------
