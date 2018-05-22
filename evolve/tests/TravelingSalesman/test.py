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
import random
from math import sqrt
import numpy as np

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.do.commandline import Command
from pts.evolve.core import reference
from pts.evolve.optimize.optimizer import show_best
from pts.evolve.genomes.list1d import G1DList
from pts.evolve.core.allele import GAlleles, GAlleleList
from pts.evolve.core.crossovers import G1DListCrossoverOX
from pts.evolve.core.mutators import G1DListMutatorSwap
from pts.core.basics.animation import Animation
from pts.evolve.tests.TravelingSalesman.plot import Plotter
from pts.core.tools.random import setup_prng, skirt_seed
from pts.core.test.implementation import TestImplementation

# -----------------------------------------------------------------

this_path = fs.absolute_path(inspect.stack()[0][1])
this_dir_path = fs.directory_of(this_path)

# -----------------------------------------------------------------

description = "solving the traveling salesman problem"

# -----------------------------------------------------------------

class TSPTest(TestImplementation):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(TSPTest, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Set properties
        self.set_properties()

        # Optimize
        self.optimize()

        # Make animation
        self.make_animation()

        # Run the reference
        self.reference()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :param:
        """

        # Call the setup function of the base class
        super(TSPTest, self).setup(**kwargs)

        # Set the random state
        setup_prng(skirt_seed)

    # -----------------------------------------------------------------

    def set_properties(self):

        """
        This function ...
        :return:
        """

        # Fix the number of cities for the TSP
        number_of_cities = 30
        width = 600
        height = 400

        # -----------------------------------------------------------------

        # Define properties
        # nparameters = 20
        nindividuals = 80
        # parameter_range = IntegerRange(0, 10)
        best_raw_score = 0.0
        # round_decimal = None
        ngenerations = 1000
        mutation_rate = 0.03
        crossover_rate = 1.0
        stats_freq = 100
        # mutation_method = "range" # or gaussian, or binary
        min_or_max = "minimize"

    # -----------------------------------------------------------------

    def optimize(self):

        """
        This function ...
        :return:
        """

        # Create the cities
        coords, cm = create_cities()

        # Create the genome
        genome = create_genome(coords)

        # Settings
        settings_optimize = dict()
        settings_optimize["output"] = None
        # settings_optimize["nparameters"] = nparameters
        settings_optimize["nindividuals"] = nindividuals
        # settings_optimize["parameter_range"] = parameter_range
        settings_optimize["best_raw_score"] = best_raw_score
        # settings_optimize["round_decimal"] = round_decimal
        settings_optimize["ngenerations"] = ngenerations
        settings_optimize["mutation_rate"] = mutation_rate
        settings_optimize["crossover_rate"] = crossover_rate
        settings_optimize["stats_freq"] = stats_freq
        # settings_optimize["mutation_method"] = mutation_method
        settings_optimize["min_or_max"] = min_or_max

        # Input
        input_optimize = dict()
        input_optimize["genome"] = genome
        input_optimize["evaluator"] = eval_func
        input_optimize["initializator"] = G1DListTSPInitializator
        input_optimize["mutator"] = G1DListMutatorSwap
        input_optimize["crossover"] = G1DListCrossoverOX

        # Create the generations plotter
        plotter = Plotter(coordinates=coords, frequency=10)

        # Set the plotter
        input_optimize["generations_plotter"] = plotter

        # Create dictionary for extra arguments to the evalutor function
        input_optimize["evaluator_kwargs"] = {"distances": cm}

        # Construct the command
        optimize = Command("optimize_continuous", "solving the traveling salesman problem", settings_optimize,
                           input_optimize, cwd=".", finish=finish_optimize)

        # Add the command
        #commands.append(optimize)

        optimizer = self.run_command(optimize)

    # -----------------------------------------------------------------

    def make_animation(self):

        """
        This function ...
        :return:
        """

        # Determine path
        #temp_path = self.optimizer.config.path

        # Determine plot path
        plot_path = fs.join(self.path, "plot")

        # Create animated gif
        animation = Animation()

        # Open the images present in the directory
        for path in fs.files_in_path(plot_path, extension="png", exact_not_name="best"):
            # Add frame to the animation
            animation.add_frame_from_file(path)

        # Save the animation
        animation_path = fs.join(self.path, "animation.gif")
        animation.saveto(animation_path)

    # -----------------------------------------------------------------

    def reference(self):

        """
        This function ...
        :return:
        """

        # Set the random state
        setup_prng(skirt_seed)

        # Make new genome, made from the original Pyevolve classes
        genome = create_genome_reference(coords)
        input_optimize["genome"] = genome

        # Solve the problem with the original Pyevolve implementation
        best = reference.call(settings_optimize, input_optimize)

        # Show the best individual
        show_best(best)

# -----------------------------------------------------------------

def cartesian_matrix(coords):

   """
   A distance matrix
   :param coords:
   """

   matrix = {}

   for i,(x1,y1) in enumerate(coords):
      for j,(x2,y2) in enumerate(coords):

         dx, dy = x1-x2,y1-y2
         dist = sqrt(dx*dx + dy*dy)
         matrix [i,j] = dist

   return matrix

# -----------------------------------------------------------------

def create_coordinates(x, y):

    """
    This function ...
    :param x:
    :param y:
    :return:
    """

    coords = []

    for i in range(len(x)):

        coords.append((x[i], y[i]))

    return coords

# -----------------------------------------------------------------

def read_coords(coord_file):

   """
   Read the coords from file
   """

   coords=[]
   for line in coord_file:

      x,y=line.strip().split(",")
      coords.append((float(x),float(y)))

   return coords

# -----------------------------------------------------------------

def tour_length(matrix, tour):

   """
   Returns the total length of the tour
   """

   total=0
   num_cities=len(tour)

   for i in range(num_cities):

      j=(i+1)%num_cities
      city_i=tour[i]
      city_j=tour[j]
      total+=matrix[city_i,city_j]

   return total

# -----------------------------------------------------------------

def G1DListTSPInitializator(genome, **kwargs):

    """
    The initializator for the TSP
    :param genome:
    """

    genome.clearList()

    list_size = genome.genomeSize
    #list_size = len(genome)

    #lst = [i for i in xrange(genome.listSize)]
    lst = [i for i in xrange(list_size)]

    #print(len(genome))
    list_size = genome.genomeSize

    #for i in xrange(genome.listSize):
    for i in xrange(list_size):

        choice = random.choice(lst)
        lst.remove(choice)
        genome.append(choice)

# -----------------------------------------------------------------

def eval_func(chromosome, **kwargs):

    """
    The evaluation function
    """

    # Get the distance matrix
    cm = kwargs.pop("distances")

    # Calculate and return the length of the city tour
    return tour_length(cm, chromosome)

# -----------------------------------------------------------------

def generate_random(ncities, xmax=800, ymax=600):

    """
    This function generates random city positions
    :param ncities:
    :param xmax:
    :param ymax:
    :return:
    """

    x = []
    y = []

    for i in xrange(ncities):

        x.append(random.randint(0, xmax))
        y.append(random.randint(0, ymax))

    # Return the x and y coordinates as numpy array
    return np.array(x), np.array(y)

# -----------------------------------------------------------------

def write_random(filename, cities, xmax=800, ymax=600):

   """
   Write random cities/positions to a text file
   """

   filehandle = open(filename, "w")
   for i in xrange(cities):
      x = random.randint(0, xmax)
      y = random.randint(0, ymax)
      filehandle.write("%d,%d\n" % (x,y))
   filehandle.close()

# -----------------------------------------------------------------

def create_cities():

    """
    This function ...
    :return:
    """

    # Generate random coordinates
    x, y = generate_random(number_of_cities, width, height)
    coords = create_coordinates(x, y)

    # Create distance matrix
    cm = cartesian_matrix(coords)

    # Return
    return coords, cm

# -----------------------------------------------------------------

def create_genome(coords):

    """
    This function ...
    :param coords
    :return:
    """

    # Set the alleles to the cities numbers
    setOfAlleles = GAlleles(homogeneous=True)
    lst = [i for i in xrange(len(coords))]
    a = GAlleleList(lst)
    setOfAlleles.add(a)

    # Create genome
    genome = G1DList(len(coords))
    genome.setParams(allele=setOfAlleles)

    # Return the genome
    return genome

# -----------------------------------------------------------------

def create_genome_reference(coords):

    """
    This function creates the genome, but using the original Pyevolve classes
    :param coords:
    :return:
    """

    from pyevolve.GAllele import GAlleles as _GAlleles
    from pyevolve.GAllele import GAlleleList as _GAlleleList
    from pyevolve.G1DList import G1DList as _G1DList

    # Set the alleles to the cities numbers
    setOfAlleles = _GAlleles(homogeneous=True)
    lst = [i for i in xrange(len(coords))]
    a = _GAlleleList(lst)
    setOfAlleles.add(a)

    # Create genome
    genome = _G1DList(len(coords))
    genome.setParams(allele=setOfAlleles)

    # Return the genome
    return genome

# -----------------------------------------------------------------
# TEST FUNCTION
# -----------------------------------------------------------------

def test(temp_path):

    """
    This function ...
    """

    pass

# -----------------------------------------------------------------
