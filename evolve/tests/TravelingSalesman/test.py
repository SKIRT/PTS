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
from functools import partial

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.do.commandline import Command
from pts.evolve import reference
from pts.evolve.optimize import show_best
from pts.evolve.genomes.list1d import G1DList
from pts.evolve.allele import GAlleles, GAlleleList
from pts.evolve.crossovers import G1DListCrossoverOX
from pts.evolve.mutators import G1DListMutatorSwap

PIL_SUPPORT = None
try:
   from PIL import Image, ImageDraw, ImageFont
   PIL_SUPPORT = True
except: PIL_SUPPORT = False

# -----------------------------------------------------------------

this_path = fs.absolute_path(inspect.stack()[0][1])
this_dir_path = fs.directory_of(this_path)

# -----------------------------------------------------------------

description = "solving the traveling salesman problem"

# -----------------------------------------------------------------

# Fix the number of cities for the TSP
number_of_cities = 30
width = 600
height = 400

# -----------------------------------------------------------------

# Define properties
#nparameters = 20
nindividuals = 80
#parameter_range = IntegerRange(0, 10)
best_raw_score = 0.0
#round_decimal = None
ngenerations = 1000
mutation_rate = 0.03
crossover_rate = 1.0
stats_freq = 100
#mutation_method = "range" # or gaussian, or binary
min_or_max = "minimize"

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

def write_tour_to_img(coords, tour, img_file):

   """
   The function to plot the graph
   :param coords:
   :param tour:
   :param img_file:
   """

   padding = 20
   coords = [(x+padding,y+padding) for (x,y) in coords]
   maxx, maxy = 0, 0

   for x, y in coords:

      maxx = max(x, maxx)
      maxy = max(y, maxy)

   maxx += padding
   maxy += padding

   img = Image.new("RGB",(int(maxx),int(maxy)),color=(255,255,255))

   font = ImageFont.load_default()
   d = ImageDraw.Draw(img)
   num_cities = len(tour)

    # Loop over the cities
   for i in range(num_cities):

      j = (i+1) % num_cities
      city_i = tour[i]
      city_j = tour[j]
      x1,y1 = coords[city_i]
      x2,y2 = coords[city_j]

      d.line((int(x1),int(y1),int(x2),int(y2)),fill=(0,0,0))
      d.text((int(x1)+7,int(y1)-5),str(i),font=font,fill=(32,32,32))

   for x, y in coords:

      x, y = int(x),int(y)
      d.ellipse((x-5,y-5,x+5,y+5),outline=(0,0,0),fill=(196,196,196))

   del d

   img.save(img_file, "PNG")

   #print("The plot was saved into the %s file." % (img_file,))

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

def evolve_callback(ga_engine, **kwargs):

    """
    This function ...
    :param ga_engine:
    :param kwargs:
    :return:
    """

    # Get the coords
    coords = kwargs.pop("coordinates")

    if ga_engine.currentGeneration % 10 == 0:
        best = ga_engine.bestIndividual()
        write_tour_to_img(coords, best, "tsp_result_%d.png" % (ga_engine.currentGeneration,))
    return False

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

# For database ...
#sqlite_adapter = DBAdapters.DBSQLite(identify="tsp", commit_freq=1000, frequency=500)

# -----------------------------------------------------------------

# Initialize list for the commands
commands = []

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

# Create the cities
coords, cm = create_cities()

# Create the genome
genome = create_genome(coords)

# Settings
settings_optimize = dict()
settings_optimize["output"] = None
#settings_optimize["nparameters"] = nparameters
settings_optimize["nindividuals"] = nindividuals
#settings_optimize["parameter_range"] = parameter_range
settings_optimize["best_raw_score"] = best_raw_score
#settings_optimize["round_decimal"] = round_decimal
settings_optimize["ngenerations"] = ngenerations
settings_optimize["mutation_rate"] = mutation_rate
settings_optimize["crossover_rate"] = crossover_rate
settings_optimize["stats_freq"] = stats_freq
#settings_optimize["mutation_method"] = mutation_method
settings_optimize["min_or_max"] = min_or_max

# Other
settings_optimize["progress_bar"] = True

# Input
input_optimize = dict()
input_optimize["genome"] = genome
input_optimize["evaluator"] = eval_func
input_optimize["initializator"] = G1DListTSPInitializator
input_optimize["mutator"] = G1DListMutatorSwap
input_optimize["crossover"] = G1DListCrossoverOX
input_optimize["callback"] = evolve_callback
#input_optimize["adapter"] = sqlite_adapter

# Create dictionary for extra arguments to the evalutor function
input_optimize["evaluator_kwargs"] = {"distances": cm}
input_optimize["callback_kwargs"] = {"coordinates": coords}

# -----------------------------------------------------------------

def finish_optimize(optimizer, **kwargs):

    """
    This function ...
    :param optimizer:
    :param kwargs:
    :return:
    """

    # Get the best
    best = optimizer.best

    # Get the coordinates
    coords = kwargs.pop("coordinates")

    # Determine path
    temp_path = optimizer.config.path
    filepath = fs.join(temp_path, "tsp_result.png")

    if PIL_SUPPORT: write_tour_to_img(coords, best, filepath)
    else: print("No PIL detected, cannot plot the graph !")

# -----------------------------------------------------------------

# Add the finish_optimize command
partialized = partial(finish_optimize, **{"coordinates": coords})
#commands.append(partialized)

# -----------------------------------------------------------------

# Construct the command
optimize = Command("optimize", "solving the traveling salesman problem", settings_optimize, input_optimize, cwd=".", finish=partialized)

# Add the command
commands.append(optimize)

# -----------------------------------------------------------------
# TEST FUNCTION
# -----------------------------------------------------------------

def test(temp_path):

    """
    This function ...
    """

    # Solve the problem with the original Pyevolve implementation
    best = reference.call(settings_optimize, input_optimize)

    # Show the best individual
    show_best(best)

# -----------------------------------------------------------------
