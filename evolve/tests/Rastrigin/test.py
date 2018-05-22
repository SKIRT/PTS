#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import math
import inspect
import numpy as np
from matplotlib import pyplot as plt

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.do.commandline import Command
from pts.core.basics.range import RealRange
from pts.evolve.core import reference
from pts.evolve.optimize.optimizer import show_best
from pts.core.test.implementation import TestImplementation

# -----------------------------------------------------------------

this_path = fs.absolute_path(inspect.stack()[0][1])
this_dir_path = fs.directory_of(this_path)

# -----------------------------------------------------------------

description = "optimizing the Rastrigin function, a deceptive function"

# -----------------------------------------------------------------

# Define properties
#nparameters = 20
nparameters = 2
nindividuals = 80
parameter_range = RealRange(-5.2, 5.30)
best_raw_score = 0.0
round_decimal = None
#ngenerations = 800
ngenerations = 100
mutation_rate = 0.05
crossover_rate = None
stats_freq = 50
mutation_method = "gaussian" # or range, or binary
min_or_max = "minimize"

# -----------------------------------------------------------------

class RastringinTest(TestImplementation):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(RastringinTest, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Optimize
        self.optimize()

        # Use pyevolve
        self.reference()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        #
        super(RastringinTest, self).setup(**kwargs)

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
        input_optimize["evaluator"] = rastringin

        # Create plot
        ax = create_plot()

        # Create callback to plot the best individuals for each generation
        input_optimize["callback"] = add_best_to_plot
        input_optimize["callback_kwargs"] = {"ax": ax}

        # Construct the command
        optimize = Command("optimize", "optimize the Rastrigin function", settings_optimize, input_optimize, cwd=".",
                           finish=show_plot)

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

def rastringin_plot(x, y):

    """
    This function ...
    :param x:
    :param y:
    :return:
    """

    n = dimension = 2

    # Calculate ...
    #result = []
    #for i in range(len(x)):
    #    xi = x[i]
    #    yi = y[i]
    #    zi = xi**2 - 10 * math.cos(2 * math.pi * xi) + yi**2 - 10 * math.cos(2 * math.pi * yi)
    #    result.append(zi)
    #result = np.array(result)

    result = x**2 - 10 * np.cos(2 * math.pi * x) + y**2 - 10 * np.cos(2 * math.pi * y)

    return (10 * dimension) + result

# -----------------------------------------------------------------

def rastringin(genome):

    """
    This is the Rastringin function, a deceptive function
    :param genome:
    :return:
    """

    n = len(genome)
    total = 0
    for i in xrange(n):
        total += genome[i] ** 2 - 10 * math.cos(2 * math.pi * genome[i])
    return (10 * n) + total

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

    #if index % 50 != 0: return

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

    numbins = 50
    numzbins = 30

    # make a contour plot of the rosenbrock function surface.
    X, Y = np.meshgrid(np.linspace(parameter_range.min, parameter_range.max, numbins), np.linspace(parameter_range.min, parameter_range.max, numbins))

    # Create Z values
    Z = rastringin_plot(X, Y)

    #ax.plot([1], [1], 'x', mew=3, markersize=10, color='#111111')

    #ax.contourf(X, Y, Z, np.logspace(-1, 3, numzbins), cmap='jet')
    ax.contourf(X, Y, Z, cmap='jet')

    # Set axes limits
    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-0.9, 1.7)

    # Return the axes
    return ax

# -----------------------------------------------------------------

def show_plot(*args, **kwargs):

    """
    This function ...
    :return:
    """

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
