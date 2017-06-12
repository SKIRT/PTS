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
from pts.core.tools import formatting as fmt
from pts.core.tools.stringify import tostr

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
definition.add_positional_optional("generations", "string_list", "name of the generations for which to show the crossover")

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

print(fitting_run.parameter_base_types)

# -----------------------------------------------------------------

# Get generation names
generations = config.generations if config.generations is not None else fitting_run.genetic_generations

# -----------------------------------------------------------------

print("")

# Loop over the generations
for generation_name in generations:

    print(generation_name)
    print("")

    # Get the generation
    generation = fitting_run.get_generation(generation_name)
    platform = fitting_run.get_generation_platform(generation_name)

    # -----------------------------------------------------------------

    # Load the crossover table
    crossover = generation.crossover_table

    # -----------------------------------------------------------------

    # Load parent and new populations
    parents = generation.parents
    newborns = generation.newborns

    # -----------------------------------------------------------------

    units = fitting_run.parameter_units

    # -----------------------------------------------------------------

    print("")
    for index in range(len(crossover)):

        mother_name, father_name = crossover.get_parents(index)
        sister_name, brother_name = crossover.get_children(index)

        mother = platform.make_genome(parents[mother_name])
        father = platform.make_genome(parents[father_name])

        sister = platform.make_genome(newborns[sister_name])
        brother = platform.make_genome(newborns[brother_name])

        mother_parameters = platform.genome_to_parameters(mother)
        father_parameters = platform.genome_to_parameters(father)
        sister_parameters = platform.genome_to_parameters(sister)
        brother_parameters = platform.genome_to_parameters(brother)

        # To strings
        mother_string = tostr(mother.genes, delimiter=" ")
        father_string = tostr(father.genes, delimiter=" ")

        # Crossover happened
        if crossover.is_crossover(index):

            # Get crossover method and details
            method = crossover.get_crossover_method(index)
            details = crossover.get_crossover_details(index)

            # Create crossover genomes
            initial_sister, initial_brother = platform.crossover(mother, father, details)

            # Check where mutation happened
            sister_mutations = [initial_sister.genes[i] != sister.genes[i] for i in range(len(initial_sister))]
            brother_mutations = [initial_brother.genes[i] != brother.genes[i] for i in range(len(initial_brother))]

            # COLORED
            mother_colored = fmt.colored_sequence(mother, colors="green", delimiter=" ")
            father_colored = fmt.colored_sequence(father, colors=None, delimiter=" ")

            # Get the origins
            sister_origins, brother_origins = platform.crossover_origins(len(mother), details)

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

            print("#" + str(index+1) + " " + fmt.blue + fmt.underlined + method.title() + " Crossover:" + fmt.reset)
            print("")
            print("Parents   :   " + mother_colored + "  x  " + father_colored)
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

            print("#" + str(index+1) + " " + fmt.blue + fmt.underlined + "Cloning:" + fmt.reset)
            print("")
            print("Parents  :    " + mother_colored + fmt.reset + "     " + father_colored)
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
