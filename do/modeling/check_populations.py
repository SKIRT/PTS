#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.evolve.check_populations Check the populations data.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.fitting.component import get_run_names, get_populations
from pts.modeling.fitting.run import get_generations_table, get_generation_path
from pts.modeling.fitting.generation import Generation
from pts.modeling.core.environment import verify_modeling_cwd

# -----------------------------------------------------------------

modeling_path = verify_modeling_cwd()

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()
run_names = get_run_names(modeling_path)
if len(run_names) == 0: raise RuntimeError("There are no fitting runs")
elif len(run_names) == 1: definition.add_fixed("fitting_run", "string", run_names[0])
else: definition.add_required("fitting_run", "string", "name of the fitting run to use", choices=run_names)

# Create the configuration
config = parse_arguments("check_populations", definition)

# -----------------------------------------------------------------

# Load populations table
populations = get_populations(modeling_path)
populations_run = populations[config.fitting_run]

# Load generation table
generations = get_generations_table(modeling_path, config.fitting_run)

# -----------------------------------------------------------------

generation_names = generations.genetic_generations_with_initial

# -----------------------------------------------------------------

newborns_dict = OrderedDict()

# -----------------------------------------------------------------

# Loop over the generations
for index in range(len(populations_run)):

    population = populations_run[index]
    generation_name = generation_names[index]

    # Get path
    generation_path = get_generation_path(modeling_path, config.fitting_run, generation_name)

    # Load generation
    generation = Generation.from_path(generation_path)
    newborns = generation.newborns
    parents = generation.parents

    #population_path = fs.join(generation_path, "population.dat")
    #newborns = load_population(population_path)

    #print(generation_name)
    #print(newborns == population)

    #for individual_key in population:
    #    individual = population[individual_key]
    #    newborn = newborns[individual_key]
    #    print(individual == newborn)

    #if newborns is not None: print(containers.equal_dictionaries(population, newborns))
    #else: print(containers.equal_dictionaries(population, parents))

    if newborns is not None: newborns_dict[generation_name] = newborns

# -----------------------------------------------------------------

for index in range(len(populations_run)):

    population = populations_run[index]
    generation_name = generation_names[index]

    if generation_name not in newborns_dict: continue

    print(generation_name)

    for individual_key in population:
        individual = population[individual_key]
        newborn = newborns_dict[generation_name][individual_key]
        print(individual == newborn) # TODO: IF NOT EQUAL: DOES THIS CORRESPOND WITH ELITISM???

# -----------------------------------------------------------------
