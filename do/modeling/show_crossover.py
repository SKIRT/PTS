#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.evolve.show_crossover Show the crossover.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.core.tools.logging import setup_log
from pts.core.tools import filesystem as fs
from pts.modeling.fitting.component import get_run_names
from pts.modeling.core.environment import GalaxyModelingEnvironment
from pts.modeling.fitting.component import load_fitting_run
from pts.evolve.optimize.optimizer import gray_binary_string_to_parameters, binary_string_to_parameters

# -----------------------------------------------------------------

modeling_path = fs.cwd()

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()
run_names = get_run_names(modeling_path)
if len(run_names) == 0: raise RuntimeError("There are no fitting runs")
elif len(run_names) == 1: definition.add_fixed("fitting_run", "string", run_names[0])
else: definition.add_required("fitting_run", "string", "name of the fitting run to use", choices=run_names)

# Generation
definition.add_required("generation_name", "string", "name of the generation for which to show the crossover")

# Create the configuration
setter = ArgumentConfigurationSetter("show_crossover")
config = setter.run(definition)

# Set logging
log = setup_log("DEBUG")

# -----------------------------------------------------------------

# Load the modeling environment
environment = GalaxyModelingEnvironment(modeling_path)

# -----------------------------------------------------------------

# Load the fitting run
fitting_run = load_fitting_run(modeling_path, config.fitting_run)

# -----------------------------------------------------------------

# Get the generation
generation = fitting_run.get_generation(config.generation_name)

# -----------------------------------------------------------------

# Load the crossover table
crossover = generation.crossover_table

# -----------------------------------------------------------------

# Load parent and new populations
parents = generation.parents
newborns = generation.newborns

# -----------------------------------------------------------------

# Load the optimizer input
optimizer_input = generation.optimizer_input

# -----------------------------------------------------------------

parameter_minima_scalar = optimizer_input["minima"]
parameter_maxima_scalar = optimizer_input["maxima"]

# -----------------------------------------------------------------

ndigits = optimizer_input["ndigits"]
nbits = optimizer_input["nbits"]
scales = optimizer_input["scales"]

# -----------------------------------------------------------------

# Load optimizer config
optimizer_config = generation.optimizer_config

# Get settings
elitism = optimizer_config.elitism
gray_code = optimizer_config.gray_code

# -----------------------------------------------------------------

units = fitting_run.parameter_units

# -----------------------------------------------------------------

print("")
for index in range(len(crossover)):

    mother_name, father_name = crossover.get_parents(index)
    sister_name, brother_name = crossover.get_children(index)

    # Get genomes of parents
    mother = parents[mother_name]
    father = parents[father_name]

    # Get genomes of children
    sister = newborns[sister_name]
    brother = newborns[brother_name]

    # Convert
    if gray_code:

        mother_real = gray_binary_string_to_parameters(mother, parameter_minima_scalar, parameter_maxima_scalar, fitting_run.nbits_list)
        father_real = gray_binary_string_to_parameters(father, parameter_minima_scalar, parameter_maxima_scalar, fitting_run.nbits_list)
        sister_real = gray_binary_string_to_parameters(sister, parameter_minima_scalar, parameter_maxima_scalar, fitting_run.nbits_list)
        brother_real = gray_binary_string_to_parameters(brother, parameter_minima_scalar, parameter_maxima_scalar, fitting_run.nbits_list)

    else:

        mother_real = binary_string_to_parameters(mother, parameter_minima_scalar, parameter_maxima_scalar, fitting_run.nbits_list)
        father_real = binary_string_to_parameters(father, parameter_minima_scalar, parameter_maxima_scalar, fitting_run.nbits_list)
        sister_real = binary_string_to_parameters(sister, parameter_minima_scalar, parameter_maxima_scalar, fitting_run.nbits_list)
        brother_real = binary_string_to_parameters(brother, parameter_minima_scalar, parameter_maxima_scalar, fitting_run.nbits_list)

    # Crossover happened
    if crossover.is_crossover(index):

        print("CROSSOVER:")
        print("")
        print("Parents:    " + str(mother) + " x " + str(father))
        print("")
        print("Children:   " + str(sister) + "   " + str(brother))
        print("")

        print("")
        print("Parents:    " + str(mother_real) + " x " + str(father_real))
        print("")
        print("Children:   " + str(sister_real) + "   " + str(brother_real))
        print("")

    # Just cloned
    else:

        print("CLONING:")
        print("")
        print("Parents:    " + str(mother) + "  " + str(father))
        print("")
        print("Children:   " + str(sister) + "  " + str(brother))
        print("")

        print("")
        print("Parents:    " + str(mother_real) + "   " + str(father_real))
        print("")
        print("Children:   " + str(sister_real) + "   " + str(brother_real))

# -----------------------------------------------------------------
