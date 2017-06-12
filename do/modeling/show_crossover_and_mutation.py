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

    units = fitting_run.parameter_units

    # -----------------------------------------------------------------

    print("")
    # Enact each reproduction event
    for reproduction in platform.reproductions:

        # Get actual parameter values for parents and children
        mother_parameters = platform.genome_to_parameters(reproduction.mother)
        father_parameters = platform.genome_to_parameters(reproduction.father)
        sister_parameters = platform.genome_to_parameters(reproduction.sister)
        brother_parameters = platform.genome_to_parameters(reproduction.brother)

        # To strings
        mother_string = tostr(reproduction.mother.genes, delimiter=" ")
        father_string = tostr(reproduction.father.genes, delimiter=" ")

        # Crossover happened
        if reproduction.crossover:

            # COLORED
            mother_colored = fmt.colored_sequence(reproduction.mother, colors="green", delimiter=" ")
            father_colored = fmt.colored_sequence(reproduction.father, colors=None, delimiter=" ")

            initial_sister_colors = ["green" if flag else None for flag in reproduction.sister_origins]
            initial_brother_colors = [None if flag else "green" for flag in reproduction.brother_origins]

            # COLORED
            initial_sister_colored = fmt.colored_sequence(reproduction.initial_sister.genes, colors=initial_sister_colors, delimiter=" ")
            initial_brother_colored = fmt.colored_sequence(reproduction.initial_brother.genes, colors=initial_brother_colors, delimiter=" ")

            # COLORED
            sister_colors = [initial_sister_colors[i] if not reproduction.sister_mutations[i] else "red" for i in range(len(reproduction.sister))]
            brother_colors = [initial_brother_colors[i] if not reproduction.brother_mutations[i] else "red" for i in range(len(reproduction.brother))]

            sister_colored = fmt.colored_sequence(reproduction.sister.genes, colors=sister_colors, delimiter=" ")
            brother_colored = fmt.colored_sequence(reproduction.brother.genes, colors=brother_colors, delimiter=" ")

            print("#" + str(reproduction.index+1) + " " + fmt.blue + fmt.underlined + generation.crossover_method.title() + " Crossover:" + fmt.reset)
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
            print("  - sister: " + str(reproduction.nmutations_sister) + " of " + str(len(reproduction.sister)) + " (" + str(reproduction.relative_nmutations_sister) + ")")
            print("  - brother: " + str(reproduction.nmutations_brother) + " of " + str(len(reproduction.brother)) + " (" + str(reproduction.relative_nmutations_brother) + ")")
            print("")

        # Just cloned
        else:

            # COLORED
            mother_colored = fmt.colored_sequence(reproduction.mother, colors="green", delimiter=" ")
            father_colored = fmt.colored_sequence(reproduction.father, colors=None, delimiter=" ")

            # COLORED
            initial_sister_colored = fmt.colored_sequence(reproduction.initial_sister, colors="green", delimiter=" ")
            initial_brother_colored = fmt.colored_sequence(reproduction.initial_brother, colors=None, delimiter=" ")

            sister_colors = ["green" if not reproduction.sister_mutations[i] else "red" for i in range(len(reproduction.sister))]
            brother_colors = [None if not reproduction.brother_mutations[i] else "red" for i in range(len(reproduction.brother))]

            sister_colored = fmt.colored_sequence(reproduction.sister.genes, colors=sister_colors, delimiter=" ")
            brother_colored = fmt.colored_sequence(reproduction.brother.genes, colors=brother_colors, delimiter=" ")

            print("#" + str(reproduction.index+1) + " " + fmt.blue + fmt.underlined + "Cloning:" + fmt.reset)
            print("")
            print("Parents  :    " + mother_colored + fmt.reset + "     " + father_colored)
            print("Cloning  :    " + initial_sister_colored + "     " + initial_brother_colored)
            print("Mutation :    " + sister_colored + "     " + brother_colored)
            print("")

            print("Number of mutations:")
            print("")
            print("  - sister: " + str(reproduction.nmutations_sister) + " of " + str(len(reproduction.sister)) + " (" + str(reproduction.relative_nmutations_sister) + ")")
            print("  - brother: " + str(reproduction.nmutations_brother) + " of " + str(len(reproduction.brother)) + " (" + str(reproduction.relative_nmutations_brother) + ")")
            print("")

            #print("")
            #print("Parents:    " + str(mother_real) + "   " + str(father_real))
            #print("")
            #print("Children:   " + str(sister_real) + "   " + str(brother_real))
            #print("")

# -----------------------------------------------------------------
