#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.platform Contains the GenerationPlatform class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from .generation import Generation
from .evaluate import get_parameter_values_from_genome
from .reproduction import ReproductionEvent
from ...core.tools import formatting as fmt
from ...core.tools.stringify import tostr

# -----------------------------------------------------------------

class GenerationPlatform(object):
    
    """
    This class...
    """

    def __init__(self, generation):

        """
        The constructor ...
        :param generation:
        :return:
        """

        self.generation = generation

    # -----------------------------------------------------------------

    @classmethod
    def from_path(cls, generation_path):

        """
        This function ...
        :param generation_path:
        :return:
        """

        # Load generation
        generation = Generation.from_path(generation_path)

        # Create and return platform
        return cls(generation)

    # -----------------------------------------------------------------

    def make_genome(self, genes):

        """
        Tihs function ...
        :param genes:
        :return:
        """

        return self.generation.genome_class(genes=genes)

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_run(self):

        """
        This function ...
        :return:
        """

        return self.generation.fitting_run

    # -----------------------------------------------------------------

    def genome_to_parameters(self, genome):

        """
        This function ...
        :param genome:
        :return:
        """

        # genome, fitting_run, minima, maxima, nbits, parameter_scales, gray=False
        return get_parameter_values_from_genome(genome, self.fitting_run, self.generation.parameter_minima_scalar, self.generation.parameter_maxima_scalar, self.generation.nbits_list, self.generation.parameter_scales, self.generation.gray_code)

    # -----------------------------------------------------------------

    def genes_to_parameters(self, genes):

        """
        This function ...
        :param genes:
        :return:
        """

        genome = self.make_genome(genes)
        return self.genome_to_parameters(genome)

    # -----------------------------------------------------------------

    def index_for_parameter(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        return self.generation.index_for_parameter(label)

    # -----------------------------------------------------------------

    def parameter_to_string(self, label, value):

        """
        This function ...
        :param label:
        :param value:
        :return:
        """

        return tostr(value, scientific=True, fancy=True, ndigits=self.generation.parameter_ndigits[self.index_for_parameter(label)])

    # -----------------------------------------------------------------

    def parameter_values_to_strings(self, values):

        """
        This function ...
        :param values: DICT
        :return:
        """

        strings = dict()

        for label in values:

            # Add
            string = self.parameter_to_string(label, values[label])
            strings[label] = string

        # Return
        return strings

    # -----------------------------------------------------------------

    def parameter_values_to_string(self, values, delimiter=", "):

        """
        This function ...
        :param values:
        :param delimiter:
        :return:
        """

        strings = self.parameter_values_to_strings(values)

        parts = []
        for label in self.generation.parameter_labels:
            parts.append(label + " = " + strings[label])

        # Make one long string
        return delimiter.join(parts)

    # -----------------------------------------------------------------

    def genome_to_parameters_string(self, genome, delimiter=", "):

        """
        This function ...
        :param genome:
        :param delimiter:
        :return:
        """

        values = self.genome_to_parameters(genome)
        return self.parameter_values_to_string(values, delimiter=delimiter)

    # -----------------------------------------------------------------

    def crossover(self, mother, father, details=None):

        """
        This function ...
        :param mother:
        :param father:
        :param details:
        :return:
        """

        sister, brother = self.generation.crossover_function(None, mom=mother, dad=father, details=details, count=2)
        return sister, brother

    # -----------------------------------------------------------------

    def crossover_origins(self, size, details):

        """
        This function ...
        :param size:
        :param details:
        :return:
        """

        sister_origins, brother_origins = self.generation.crossover_origins_function(size, details)
        return sister_origins, brother_origins

    # -----------------------------------------------------------------

    @property
    def nreproductions(self):

        """
        This function ...
        :return:
        """

        return len(self.generation.crossover_table)

    # -----------------------------------------------------------------

    @property
    def ncrossovers(self):

        """
        This function ...
        :return:
        """

        return self.generation.crossover_table.ncrossovers

    # -----------------------------------------------------------------

    @property
    def reproductions(self):

        """
        This function ...
        :return:
        """

        for index in range(self.nreproductions): yield self.get_reproduction(index)

    # -----------------------------------------------------------------

    def get_reproduction(self, index):

        """
        This fucntion ...
        :param index:
        :return:
        """

        # Get names of the individuals
        mother_name, father_name = self.generation.crossover_table.get_parents(index)
        sister_name, brother_name = self.generation.crossover_table.get_children(index)

        # Load parent and new populations
        parents = self.generation.parents
        newborns = self.generation.newborns

        # Create genomes of parents
        mother = self.make_genome(parents[mother_name])
        father = self.make_genome(parents[father_name])

        # Create children genomes
        sister = self.make_genome(newborns[sister_name])
        brother = self.make_genome(newborns[brother_name])

        crossover = False

        # Crossover happened
        if self.generation.crossover_table.is_crossover(index):

            crossover = True

            # Get crossover method and details
            #method = self.generation.crossover_table.get_crossover_method(index)
            details = self.generation.crossover_table.get_crossover_details(index)

            # Create crossover genomes
            initial_sister, initial_brother = self.crossover(mother, father, details)

            # Get the origins
            sister_origins, brother_origins = self.crossover_origins(len(mother), details)

        # Just cloned
        else: initial_sister, initial_brother, sister_origins, brother_origins = mother, father, None, None

        # Return
        return ReproductionEvent(index, mother, father, initial_sister, initial_brother, sister, brother, crossover, sister_origins, brother_origins)

    # -----------------------------------------------------------------

    def show_reproductions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the reproduction events ...")

        print("")

        # Enact each reproduction event
        for reproduction in self.reproductions:

            # COLORED
            mother_colored = fmt.colored_sequence(reproduction.mother, colors="green", delimiter=" ")
            father_colored = fmt.colored_sequence(reproduction.father, colors=None, delimiter=" ")

            # Crossover happened
            if reproduction.crossover:

                reproduction_type_text = "crossover"

                # Colors for initial siblings
                initial_sister_colors = ["green" if flag else None for flag in reproduction.sister_origins]
                initial_brother_colors = [None if flag else "green" for flag in reproduction.brother_origins]

                # COLORED
                initial_sister_colored = fmt.colored_sequence(reproduction.initial_sister.genes, colors=initial_sister_colors, delimiter=" ")
                initial_brother_colored = fmt.colored_sequence(reproduction.initial_brother.genes, colors=initial_brother_colors, delimiter=" ")

                # COLORED
                sister_colors = [initial_sister_colors[i] if not reproduction.sister_mutations[i] else "red" for i in range(len(reproduction.sister))]
                brother_colors = [initial_brother_colors[i] if not reproduction.brother_mutations[i] else "red" for i in range(len(reproduction.brother))]

                print("#" + str(reproduction.index+1) + " " + fmt.blue + fmt.underlined + self.generation.crossover_method.title() + " Crossover:" + fmt.reset)
                print("")

            # Just cloned
            else:

                reproduction_type_text = "cloning"

                # COLORED
                initial_sister_colored = fmt.colored_sequence(reproduction.initial_sister.genes, colors="green", delimiter=" ")
                initial_brother_colored = fmt.colored_sequence(reproduction.initial_brother.genes, colors=None, delimiter=" ")

                sister_colors = ["green" if not reproduction.sister_mutations[i] else "red" for i in range(len(reproduction.sister))]
                brother_colors = [None if not reproduction.brother_mutations[i] else "red" for i in range(len(reproduction.brother))]

                print("#" + str(reproduction.index+1) + " " + fmt.blue + fmt.underlined + "Cloning:" + fmt.reset)
                print("")

            sister_colored = fmt.colored_sequence(reproduction.sister.genes, colors=sister_colors, delimiter=" ")
            brother_colored = fmt.colored_sequence(reproduction.brother.genes, colors=brother_colors, delimiter=" ")

            print("GENOMES:")
            print("")

            # Print in columns
            with fmt.print_in_columns(5) as print_row:

                print_row("Parents", ":", mother_colored, "     ", father_colored)
                print_row(reproduction_type_text.title(), ":", initial_sister_colored, "     ", initial_brother_colored)
                print_row("Mutation", ":", sister_colored, "     ", brother_colored)

            print("")
            print("Number of mutations:")
            with fmt.itemize() as print_item:
                print_item("sister: " + str(reproduction.nmutations_sister) + " of " + str(len(reproduction.sister)) + " (" + str(reproduction.relative_nmutations_sister) + ")")
                print_item("brother: " + str(reproduction.nmutations_brother) + " of " + str(len(reproduction.brother)) + " (" + str(reproduction.relative_nmutations_brother) + ")")

            # Get actual parameter values for parents and children
            mother_string = self.genome_to_parameters_string(reproduction.mother)
            father_string = self.genome_to_parameters_string(reproduction.father)
            initial_sister_string = self.genome_to_parameters_string(reproduction.initial_sister)
            initial_brother_string = self.genome_to_parameters_string(reproduction.initial_brother)
            sister_string = self.genome_to_parameters_string(reproduction.sister)
            brother_string = self.genome_to_parameters_string(reproduction.brother)

            print("PARAMETERS:")
            print("")

            with fmt.print_in_columns(5) as print_row:

                print_row("Parents", ":", "[" + mother_string + "]", "     ", "[" + father_string + "]")
                print_row(reproduction_type_text.title(), ":", "[" + initial_sister_string + "]", "     ", "[" + initial_brother_string + "]")
                print_row("Mutation", ":", "[" + sister_string + "]", "     ", "[" + brother_string + "]")

            print("")

    # -----------------------------------------------------------------

    @property
    def nrecurrences(self):

        """
        This function ...
        :return:
        """

        return len(self.generation.recurrent_data)

    # -----------------------------------------------------------------

    @property
    def recurrences(self):

        """
        This function ...
        :return:
        """

        for index in range(self.nrecurrences): yield self.get_recurrence(index)

    # -----------------------------------------------------------------

    def get_recurrence(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        # Get the key
        individual_key = self.generation.recurrent_data.keys()[index]

        # Get the chi squared
        chi_squared = self.generation.recurrent_data[individual_key]

        return individual_key, chi_squared

    # -----------------------------------------------------------------

    def show_recurrence(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the recurrency ...")

        print("")

        # Loop over each recurrence event
        for recurrence in self.recurrences:

            print(recurrence)

    # -----------------------------------------------------------------

    @property
    def nelitisms(self):

        """
        This function ...
        :return:
        """

        return len(self.generation.elitism_table)

    # -----------------------------------------------------------------

    @property
    def elitisms(self):

        """
        This function ...
        :return:
        """

        for index in range(self.nelitisms): yield self.get_elitism(index)

    # -----------------------------------------------------------------

    def get_elitism(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def show_elitism(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the elitism events ...")

        print("")

        # Loop over each elitism event
        for elitism in self.elitisms: pass

# -----------------------------------------------------------------
