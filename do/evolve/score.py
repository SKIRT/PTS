#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.evolve.score Set scores for GA test.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.core.tools import tables

# -----------------------------------------------------------------

x = np.linspace(12,25,100)

test_data_x = [20., 16., 19.79999924, 18.39999962, 17.10000038, 15.5, 14.69999981, 17.10000038, 15.39999962,
               16.20000076,
               15., 17.20000076, 16., 17., 14.39999962]
test_data_y = [88.59999847, 71.59999847, 93.30000305, 84.30000305, 80.59999847, 75.19999695, 69.69999695, 82.,
               69.40000153, 83.30000305, 79.59999847, 82.59999847, 80.59999847, 83.5, 76.30000305]

# -----------------------------------------------------------------

def fit_function(x, a, b):

    """
    This function ...
    :param x:
    :param a:
    :param b:
    :return:
    """

    return a * x + b

# -----------------------------------------------------------------

def chi_squared_function(chromosome):

    """
    This function calculates the chi-squared value for a certain set of parameters (chromosome)
    :param chromosome:
    :return:
    """

    chi_squared = 0.0
    for i in range(len(test_data_x)):
        x = test_data_x[i]
        y = test_data_y[i]
        chromosome_y = fit_function(x, chromosome[0], chromosome[1])
        chi_squared += (y - chromosome_y) ** 2.
    chi_squared /= 2.0
    return chi_squared

# -----------------------------------------------------------------

last_generation = None

# Check the index of the last generation
for name in fs.directories_in_path():
    if "ref" in name: continue
    generation = int(name.split("Generation ")[1])
    if last_generation is None or generation > last_generation: last_generation = generation

# -----------------------------------------------------------------

if last_generation is None:
    generation_path = fs.cwd()
    print("current generation: the initial population")
else:
    generation_path = fs.join(fs.cwd(), "Generation " + str(last_generation))
    print("Current generation: " + str(last_generation))

# -----------------------------------------------------------------

# Path to the current GA object
path = fs.join(generation_path, "ga.pickle")

# Path to the parameters table
parameters_path = fs.join(generation_path, "parameters.dat")

# -----------------------------------------------------------------

# Load the parameters table
table = tables.from_file(parameters_path, format="ascii.ecsv")

# -----------------------------------------------------------------

names = []
scores = []

lowest_score = None
index_lowest = None

for index in range(len(table)):

    # Get the parameter values
    parameter_a_tab = table["Parameter a"][index]
    parameter_b_tab = table["Parameter b"][index]

    # Calculate the score
    score = chi_squared_function([parameter_a_tab, parameter_b_tab])

    # Keep track of index of lowest score
    if lowest_score is None or score < lowest_score:
        lowest_score = score
        index_lowest = index

    # Add the score to the list
    name = table["Unique name"][index]
    names.append(name)
    scores.append(score)

# Create the chi squared table
data = [names, scores]
names = ["Unique name", "Chi-squared"]
chi_squared_table = tables.new(data, names)

# Determine the path to the chi squared table
chi_squared_path = fs.join(generation_path, "chi_squared.dat")

# Write the chi squared table
tables.write(chi_squared_table, chi_squared_path, format="ascii.ecsv")

# -----------------------------------------------------------------

best_parameter_a = table["Parameter a"][index_lowest]
best_parameter_b = table["Parameter b"][index_lowest]

best_path = fs.join(generation_path, "best.dat")

with open(best_path, 'w') as best_file:
    best_file.write("Parameter a: " + str(best_parameter_a) + "\n")
    best_file.write("Parameter b: " + str(best_parameter_b) + "\n")

popt, pcov = curve_fit(fit_function, test_data_x, test_data_y)

parameter_a_real = popt[0]
parameter_b_real = popt[1]

print("Best parameter a:", best_parameter_a, " REAL:", parameter_a_real)
print("Best parameter b:", best_parameter_b, " REAL:", parameter_b_real)

plt.figure()
plt.scatter(test_data_x, test_data_y)
plt.plot(x, [fit_function(best_parameter_a, best_parameter_b, x_i) for x_i in x])
plt.ylim(65, 95)
plt.xlim(12,22)

# Save the figure
plot_path = fs.join(generation_path, "best.pdf")
plt.savefig(plot_path)

# -----------------------------------------------------------------
