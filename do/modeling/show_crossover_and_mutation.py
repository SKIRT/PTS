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
from pts.core.tools import formatting as fmt
from pts.core.tools.stringify import tostr
from pts.evolve.optimize.optimizer import list_crossovers_1d, binary_string_crossovers_1d, genomes_1d, genomes_2d
from pts.evolve.optimize.optimizer import list_crossover_origins_1d, list_crossover_origins_2d, binary_string_crossover_origins_1d, binary_string_crossover_origins_2d

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
genome_type = optimizer_config.genome_type

# -----------------------------------------------------------------

units = fitting_run.parameter_units

# -----------------------------------------------------------------

print("")
for index in range(len(crossover)):

    mother_name, father_name = crossover.get_parents(index)
    sister_name, brother_name = crossover.get_children(index)

    # Genome
    make_genome = genomes_1d[genome_type]

    # Get genomes of parents
    mother = make_genome(genes=parents[mother_name])
    father = make_genome(genes=parents[father_name])

    # Get genomes of children
    sister = make_genome(genes=newborns[sister_name])
    brother = make_genome(genes=newborns[brother_name])

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

    # To strings
    mother_string = tostr(mother.genes, delimiter=" ")
    father_string = tostr(father.genes, delimiter=" ")

    # Crossover happened
    if crossover.is_crossover(index):

        # Get crossover method and details
        method = crossover.get_crossover_method(index)
        details = crossover.get_crossover_details(index)

        # Get the correct crossover function
        if genome_type == "list": crossover_function = list_crossovers_1d[method]
        elif genome_type == "binary_string": crossover_function = binary_string_crossovers_1d[method]
        else: raise ValueError("Invalid genome type")

        # Create crossover genomes
        initial_sister, initial_brother = crossover_function(None, mom=mother, dad=father, details=details, count=2)

        # Check where mutation happened
        sister_mutations = [initial_sister.genes[i] != sister.genes[i] for i in range(len(initial_sister))]
        brother_mutations = [initial_brother.genes[i] != brother.genes[i] for i in range(len(initial_brother))]

        # Get the correct crossover origins function
        if genome_type == "list": origins_function = list_crossover_origins_1d[method]
        elif genome_type == "binary_string": origins_function = binary_string_crossover_origins_1d[method]
        else: raise ValueError("Invalid genome type")

        # COLORED
        mother_colored = fmt.colored_sequence(mother, colors="green", delimiter=" ")
        father_colored = fmt.colored_sequence(father, colors=None, delimiter=" ")

        # Get the origins
        sister_origins, brother_origins = origins_function(len(mother), details)

        initial_sister_colors = ["green" if flag else None for flag in sister_origins]
        initial_brother_colors = [None if flag else "green" for flag in brother_origins]

        # COLORED
        initial_sister_colored = fmt.colored_sequence(initial_sister.genes, colors=initial_sister_colors, delimiter=" ")
        initial_brother_colored = fmt.colored_sequence(initial_brother.genes, colors=initial_brother_colors, delimiter=" ")

        # COLORED
        sister_colors = [initial_sister_colors[i] if not sister_mutations[i] else "red" for i in range(len(sister))]
        brother_colors = [initial_brother_colors[i] if not brother_mutations[i] else "red" for i in range(len(brother))]

        sister_colored = fmt.colored_sequence(sister.genes, colors=sister_colors, delimiter=" ")
        brother_colored = fmt.colored_sequence(brother.genes, colors=brother_colors, delimiter=" ")

        nmutations_sister = sum(sister_mutations)
        nmutations_brother = sum(brother_mutations)
        relative_nmutations_sister = float(nmutations_sister) / len(sister)
        relative_nmutations_brother = float(nmutations_brother) / len(brother)

        #print(nmutations_sister, nmutations_brother)

        print(fmt.blue + fmt.underlined + method.title() + " Crossover:" + fmt.reset)
        print("")
        print("Parents   :   " + mother_colored + "  x  " + father_colored)
        #print("")
        print("Crossover :   " + initial_sister_colored + "     " + initial_brother_colored)
        print("Mutation  :   " + sister_colored + "     " + brother_colored)
        print("")

        #print("")
        #print("Parents:    " + str(mother_real) + " x " + str(father_real))
        #print("")
        #print("Children:   " + str(sister_real) + "   " + str(brother_real))
        #print("")

        print("Number of mutations:")
        print("")
        print("  - sister: " + str(nmutations_sister) + " of " + str(len(sister)) + " (" + str(relative_nmutations_sister) + ")")
        print("  - brother: " + str(nmutations_brother) + " of " + str(len(brother)) + " (" + str(relative_nmutations_brother) + ")")
        print("")

    # Just cloned
    else:

        # COLORED
        mother_colored = fmt.colored_sequence(mother, colors="green", delimiter=" ")
        father_colored = fmt.colored_sequence(father, colors=None, delimiter=" ")

        # COLORED
        initial_sister_colored = fmt.colored_sequence(mother, colors="green", delimiter=" ")
        initial_brother_colored = fmt.colored_sequence(father, colors=None, delimiter=" ")

        # Check where mutation happened
        sister_mutations = [sister.genes[i] != mother.genes[i] for i in range(len(sister))]
        brother_mutations = [brother.genes[i] != father.genes[i] for i in range(len(brother))]

        sister_colors = ["green" if not sister_mutations[i] else "red" for i in range(len(sister))]
        brother_colors = [None if not brother_mutations[i] else "red" for i in range(len(brother))]

        sister_colored = fmt.colored_sequence(sister.genes, colors=sister_colors, delimiter=" ")
        brother_colored = fmt.colored_sequence(brother.genes, colors=brother_colors, delimiter=" ")

        nmutations_sister = sum(sister_mutations)
        nmutations_brother = sum(brother_mutations)
        relative_nmutations_sister = float(nmutations_sister) / len(sister)
        relative_nmutations_brother = float(nmutations_brother) / len(brother)

        print(fmt.blue + fmt.underlined + "Cloning:" + fmt.reset)
        print("")
        print("Parents  :    " + mother_colored + fmt.reset + "     " + father_colored)
        #print("")
        print("Cloning  :    " + initial_sister_colored + "     " + initial_brother_colored)
        print("Mutation :    " + sister_colored + "     " + brother_colored)
        print("")

        print("Number of mutations:")
        print("")
        print("  - sister: " + str(nmutations_sister) + " of " + str(len(sister)) + " (" + str(relative_nmutations_sister) + ")")
        print("  - brother: " + str(nmutations_brother) + " of " + str(len(brother)) + " (" + str(relative_nmutations_brother) + ")")
        print("")

        #print("")
        #print("Parents:    " + str(mother_real) + "   " + str(father_real))
        #print("")
        #print("Children:   " + str(sister_real) + "   " + str(brother_real))
        #print("")

# -----------------------------------------------------------------
