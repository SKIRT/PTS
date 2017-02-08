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

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.do.commandline import Command
from pts.core.basics.range import IntegerRange
from pts.evolve import reference
from pts.evolve.optimize import show_best

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

# Define properties
#nparameters = 20
nindividuals = 80
#parameter_range = IntegerRange(0, 10)
#best_raw_score = 0.0
#round_decimal = None
ngenerations = 1000
mutation_rate = 0.03
crossover_rate = 1.0
stats_freq = 100
#mutation_method = "range" # or gaussian, or binary
min_or_max = "minimize"

# -----------------------------------------------------------------

def cartesian_matrix(coords):

   """ A distance matrix """

   matrix={}

   for i,(x1,y1) in enumerate(coords):
      for j,(x2,y2) in enumerate(coords):

         dx,dy=x1-x2,y1-y2
         dist=sqrt(dx*dx + dy*dy)
         matrix[i,j]=dist

   return matrix

# -----------------------------------------------------------------

def read_coords(coord_file):

   """ Read the coords from file """

   coords=[]
   for line in coord_file:

      x,y=line.strip().split(",")
      coords.append((float(x),float(y)))

   return coords

# -----------------------------------------------------------------

def tour_length(matrix, tour):

   """ Returns the total length of the tour """

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

   """ The function to plot the graph """

   padding=20
   coords=[(x+padding,y+padding) for (x,y) in coords]
   maxx,maxy=0,0
   for x,y in coords:
      maxx=max(x,maxx)
      maxy=max(y,maxy)
   maxx+=padding
   maxy+=padding

   img = Image.new("RGB",(int(maxx),int(maxy)),color=(255,255,255))

   font=ImageFont.load_default()
   d=ImageDraw.Draw(img)
   num_cities=len(tour)

   for i in range(num_cities):

      j = (i+1)%num_cities
      city_i = tour[i]
      city_j = tour[j]
      x1,y1 = coords[city_i]
      x2,y2 = coords[city_j]
      d.line((int(x1),int(y1),int(x2),int(y2)),fill=(0,0,0))
      d.text((int(x1)+7,int(y1)-5),str(i),font=font,fill=(32,32,32))

   for x,y in coords:

      x,y=int(x),int(y)
      d.ellipse((x-5,y-5,x+5,y+5),outline=(0,0,0),fill=(196,196,196))

   del d

   img.save(img_file, "PNG")

   print("The plot was saved into the %s file." % (img_file,))

# -----------------------------------------------------------------

def G1DListTSPInitializator(genome, **args):

   """ The initializator for the TSP """

   genome.clearList()
   lst = [i for i in xrange(genome.listSize)]

   for i in xrange(genome.listSize):
      choice = random.choice(lst)
      lst.remove(choice)
      genome.append(choice)

cm = []
coords = []

# -----------------------------------------------------------------

def eval_func(chromosome):

   """ The evaluation function """

   global cm
   return tour_length(cm, chromosome)

# -----------------------------------------------------------------

def write_random(filename, cities, xmax=800, ymax=600):

   """ Write random cities/positions to a text file """

   filehandle = open(filename, "w")
   for i in xrange(cities):
      x = random.randint(0, xmax)
      y = random.randint(0, ymax)
      filehandle.write("%d,%d\n" % (x,y))
   filehandle.close()

# -----------------------------------------------------------------

from pyevolve import G1DList
from pyevolve import GSimpleGA
from pyevolve import GAllele
from pyevolve import Mutators
from pyevolve import Crossovers
from pyevolve import Consts

def main_run():

    global cm, coords

    # write_random(filename, number of the cities, max width, max_height)
    write_random("tsp_coords.txt", 30, 600, 400)

    # load the tsp data file
    filehandle = open("tsp_coords.txt", "rw")
    coords = read_coords(filehandle)
    cm = cartesian_matrix(coords)

    # Set the alleles to the cities numbers
    setOfAlleles = GAllele.GAlleles(homogeneous=True)
    lst = [i for i in xrange(len(coords))]
    a = GAllele.GAlleleList(lst)
    setOfAlleles.add(a)

    # Create genome
    genome = G1DList.G1DList(len(coords))
    genome.setParams(allele=setOfAlleles)

    genome.evaluator.set(eval_func)

    genome.mutator.set(Mutators.G1DListMutatorSwap)

    genome.crossover.set(Crossovers.G1DListCrossoverOX)

    genome.initializator.set(G1DListTSPInitializator)

    ga = GSimpleGA.GSimpleGA(genome)
    #ga.setGenerations(1000)
    #ga.setMinimax(Consts.minimaxType["minimize"])
    #ga.setCrossoverRate(1.0)
    #ga.setMutationRate(0.03)
    #ga.setPopulationSize(80)

    # sqlite_adapter = DBAdapters.DBSQLite(identify="tsp", commit_freq=1000, frequency=500)
    # ga.setDBAdapter(sqlite_adapter)

    # This is to make a video
    # ga.stepCallback.set(evolve_callback)

    ga.evolve(freq_stats=100)
    best = ga.bestIndividual()
    print(best)

    if PIL_SUPPORT: write_tour_to_img(coords, best, "tsp_result.png")
    else: print("No PIL detected, cannot plot the graph !")

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

# Other
settings_optimize["progress_bar"] = True

# Input
input_optimize = dict()
input_optimize["evaluator"] = rosenbrock

# Construct the command
optimize = Command("optimize", "optimize the Rosenbrock function", settings_optimize, input_optimize, cwd=".")

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
    best = reference.call(settings_optimize, rosenbrock)

    # Show the best individual
    show_best(best)

# -----------------------------------------------------------------
