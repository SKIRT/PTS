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

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from .generation import Generation
from .evaluate import get_parameter_values_from_genome
from .reproduction import ReproductionEvent
from .recurrence import Recurrence
from .elitism import Elitism
from ...core.tools import formatting as fmt
from ...core.tools.stringify import tostr
from ...core.tools import filesystem as fs
from ...evolve.analyse.database import get_scores_named_individuals, get_individual, get_individuals_scores
from ...core.tools import sequences
from pts.core.tools.utils import lazyproperty

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

    @property
    def generation_name(self):

        """
        This function ...
        :return:
        """

        return self.generation.name

    # -----------------------------------------------------------------

    @property
    def generation_index(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.get_genetic_generation_index(self.generation_name)

    # -----------------------------------------------------------------

    @property
    def generation_index_database(self):

        """
        This function ...
        :return:
        """

        return self.generation_index + 1

    # -----------------------------------------------------------------

    @property
    def previous_generation_name(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.get_previous_genetic_generation_name(self.generation_name)

    # -----------------------------------------------------------------

    @property
    def previous_generation_index(self):

        """
        This function ...
        :return:
        """

        return self.generation_index - 1

    # -----------------------------------------------------------------

    @property
    def previous_generation_index_database(self):

        """
        This function ...
        :return:
        """

        return self.previous_generation_index + 1

    # -----------------------------------------------------------------

    @property
    def next_generation_name(self):

        """
        This function ...
        :return: 
        """

        return self.fitting_run.get_next_genetic_generation_name(self.generation_name)

    # -----------------------------------------------------------------

    @property
    def next_generation_index(self):

        """
        This function ...
        :return:
        """

        return self.generation_index + 1

    # -----------------------------------------------------------------

    @property
    def next_generation_index(self):

        """
        This function ...
        :return:
        """

        return self.next_generation_index + 1

    # -----------------------------------------------------------------

    @property
    def modeling_path(self):

        """
        This function ...
        :return:
        """

        return self.generation.modeling_path

    # -----------------------------------------------------------------

    @lazyproperty
    def population(self):

        """
        This function ...
        :return:
        """

        #from .component import get_populations

        # Load the populations data
        #populations = get_populations(self.modeling_path)
        populations = self.fitting_context.populations
        return populations[self.generation.fitting_run_name][self.generation.name]

    # -----------------------------------------------------------------

    @lazyproperty
    def database(self):

        """
        This function ...
        :return:
        """

        #from .component import get_database
        #return get_database(self.modeling_path)

        return self.fitting_context.database

    # -----------------------------------------------------------------

    @lazyproperty
    def statistics(self):

        """
        This function ...
        :return:
        """

        #from .component import get_statistics
        #return get_statistics(self.modeling_path)

        return self.fitting_context.statistics

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

    @lazyproperty
    def fitting_context(self):

        """
        This function ...
        :return:
        """

        from .context import FittingContext
        return FittingContext.from_modeling_path(self.modeling_path)

    # -----------------------------------------------------------------

    @property
    def fitting_run_name(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.name

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

    def ndigits_for_parameter(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        return self.generation.parameter_ndigits[self.index_for_parameter(label)]

    # -----------------------------------------------------------------

    def add_unit(self, value, label):

        """
        This function ...
        :param value:
        :param label:
        :return:
        """

        if not self.has_unit(label): return value
        else: return value * self.unit_for_parameter(label)

    # -----------------------------------------------------------------

    def parameter_to_string(self, label, value):

        """
        This function ...
        :param label:
        :param value:
        :return:
        """

        return tostr(value, scientific=True, fancy=True, ndigits=self.ndigits_for_parameter(label))

    # -----------------------------------------------------------------

    def unit_for_parameter(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        return self.generation.unit_for_parameter(label)

    # -----------------------------------------------------------------

    def has_unit(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        return self.unit_for_parameter(label) is not None

    # -----------------------------------------------------------------

    @property
    def parameter_units(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.parameter_units

    # -----------------------------------------------------------------

    @property
    def parameter_labels(self):

        """
        This function ...
        :return:
        """

        return self.generation.parameter_labels

    # -----------------------------------------------------------------

    def parameter_list_to_values_dict(self, parameters):

        """
        This function ...
        :param parameters:
        :return:
        """

        values = dict()
        for index, label in enumerate(self.parameter_labels):

            # Add the unit
            if self.has_unit(label) is not None: value = parameters[index] * self.unit_for_parameter(label)
            else: value = parameters[index]

            # Set the value
            values[label] = value

        # Return the dictionary
        return values

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

    def parameter_list_to_string(self, parameters, delimiter=", "):

        """
        This function ...
        :param parameters:
        :param delimiter:
        :return:
        """

        values = self.parameter_list_to_values_dict(parameters)
        return self.parameter_values_to_string(values, delimiter=delimiter)

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
        for label in self.parameter_labels:
            parts.append(label + " = " + strings[label])

        # Make one long string
        return delimiter.join(parts)

    # -----------------------------------------------------------------

    def parameter_values_to_array(self, values):

        """
        This function ...
        :param values:
        :return:
        """

        scalar = []
        for label in self.parameter_labels:
            scalar.append(values[label].to(self.generation.unit_for_parameter(label)).value)
        return np.array(scalar)

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

    def genome_to_parameters_array(self, genome):

        """
        This function ...
        :param genome:
        :return:
        """

        values = self.genome_to_parameters(genome)
        return self.parameter_values_to_array(values)

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

        return self.generation.recurrence_table.nrecurrences

    # -----------------------------------------------------------------

    @property
    def recurrences(self):

        """
        This function ...
        :return:
        """

        for index in range(self.nrecurrences): yield self.get_recurrence(index)

    # -----------------------------------------------------------------

    @lazyproperty
    def recurrence_names(self):

        """
        This function ...
        :return:
        """

        return self.generation.recurrence_table.individual_ids

    # -----------------------------------------------------------------

    def get_recurrence_score(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.generation.recurrence_table.get_score_for_individual(name)

    # -----------------------------------------------------------------

    def get_recurrence_parameters(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        parameters = self.generation.recurrence_table.get_individual(name)
        original_parameters = self.generation.recurrence_table.get_original_individual(name)
        #print(parameters, original_parameters)
        return parameters, original_parameters

    # -----------------------------------------------------------------

    def get_recurrence_generation_and_name(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        generation_index, individual_name = self.generation.recurrence_table.get_original_generation_and_individual_id(name)
        return generation_index - 1, individual_name

    # -----------------------------------------------------------------

    def get_recurrence(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        # Get the index'th individual name
        name = self.recurrence_names[index]

        # Get the generation index and the original name
        generation_index, original_name = self.get_recurrence_generation_and_name(name)

        # Get the parameters and the original parameters
        parameters, original_parameters = self.get_recurrence_parameters(name)

        # Get the generation name
        generation_name = self.fitting_run.get_genetic_generation_name(generation_index)

        # Get newborns
        newborns = self.generation.newborns

        # Create genome
        genome = self.make_genome(newborns[name])

        # Get score
        score = self.get_recurrence_score(name)

        # Get genome of original individual
        population = self.fitting_run.get_population_for_generation(generation_name)
        original_genome = self.make_genome(population[original_name])

        #print(genome, type(genome))
        #print(original_genome, type(original_genome))

        # Create and return the recurrence object
        return Recurrence(index, genome, generation_name, original_genome, score, parameters, original_parameters, name, original_name)

    # -----------------------------------------------------------------

    def show_recurrence(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the recurrence ...")

        print("")

        # Loop over each recurrence event
        for recurrence in self.recurrences:

            # COLORED
            #individual_colored = fmt.colored_sequence(recurrence.individual, colors="green", delimiter=" ")
            #original_colored = fmt.colored_sequence(recurrence.original, colors=None, delimiter=" ")

            # COLORED
            colors = [None if not recurrence.differences[i] else "magenta" for i in range(recurrence.genome_size)]
            individual_colored = fmt.colored_sequence(recurrence.individual, None, background_colors=colors, delimiter=" ")
            original_colored = fmt.colored_sequence(recurrence.original, None, background_colors=colors, delimiter=" ")

            print("#" + str(recurrence.index + 1) + " " + fmt.blue + fmt.underlined + "Recurrence:" + fmt.reset)
            print("")

            print("NAMES:")
            print("")

            with fmt.print_in_columns(4) as print_row:

                print_row("Individual", ":", recurrence.name, "[" + self.generation_name + "]")
                print_row("Original", ":", recurrence.original_name, "[" + recurrence.generation + "]")

            print("")
            print("GENOMES:")
            print("")

            with fmt.print_in_columns(3) as print_row:

                print_row("Individual", ":", individual_colored)
                print_row("Original", ":", original_colored)

            print("")
            print(" - number of differences:" + tostr(recurrence.ndifferences))

            print("")
            print("COMPARISON:")
            print("")

            parameters = recurrence.parameters
            original_parameters = recurrence.original_parameters

            # To strings
            #par_string = self.parameter_values_to_string(parameters)
            #original_par_string = self.parameter_values_to_string(original_parameters)
            par_string = self.parameter_list_to_string(parameters)
            original_par_string = self.parameter_list_to_string(original_parameters)

            with fmt.print_in_columns(3) as print_row:

                print_row("Individual", ":", "[" + par_string + "]")
                print_row("Original", ":", "[" + original_par_string + "]")

            #print("")

            print("")
            print("PARAMETERS:")
            print("")

            # Get actual parameter values for parents and children
            parameters_string = self.genome_to_parameters_string(recurrence.individual)
            original_parameters_string = self.genome_to_parameters_string(recurrence.original) # DO WE NEED TO TAKE INTO ACCOUNT THE FACT THAT THIS GENOME COMES FROM ANOTHER GENERATION FOR THE CONVERSION?

            # Print parameters
            with fmt.print_in_columns(3) as print_row:

                print_row("Individual", ":", "[" + parameters_string + "]")
                print_row("Original", ":", "[" + original_parameters_string + "]")

            print("")

            parameters_array = self.genome_to_parameters_array(recurrence.individual)
            original_parameters_array = self.genome_to_parameters_array(recurrence.original)
            ratio = parameters_array / original_parameters_array
            absolute = abs(parameters_array - original_parameters_array)
            relative = absolute / original_parameters_array

            # Add unit for absolute al values
            absolute = [self.add_unit(value, label) for value, label in zip(absolute, self.parameter_labels)]
            absolute_strings = [self.parameter_to_string(label, value) for value, label in zip(absolute, self.parameter_labels)]

            print("Difference:")
            with fmt.itemize() as print_item:
                print_item("absolute: " + tostr(absolute_strings))
                print_item("relative: " + tostr(relative))
                print_item("ratio: " + tostr(ratio))

            # Check if the values are close enough
            if not np.isclose(ratio, 1., rtol=5e-3): log.warning("Values are not close")

            #print("")

    # -----------------------------------------------------------------

    @property
    def nelitism_attempts(self):

        """
        This function ...
        :return:
        """

        return len(self.generation.elitism_table)

    # -----------------------------------------------------------------

    @property
    def nelitisms(self):

        """
        This function ...
        :return:
        """

        return self.generation.elitism_table.nelitisms

    # -----------------------------------------------------------------

    @property
    def elitisms(self):

        """
        This function ...
        :return:
        """

        for index in range(self.nelitisms): yield self.get_elitism(index)

    # -----------------------------------------------------------------

    @lazyproperty
    def elitism_names(self):

        """
        This function ...
        :return:
        """

        return self.generation.elitism_table.elitism_individual_ids

    # -----------------------------------------------------------------

    def get_elitism(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        # Get the index'th individual name
        name = self.elitism_names[index]

        # Get the score of the individual
        score = self.generation.elitism_table.get_score_for_individual(name)

        # Get the name of the replacement individual
        replacement_name = self.generation.elitism_table.get_replacement_for_individual(name)

        # Get the score of the replacement individual
        replacement_score = self.generation.elitism_table.get_replacement_score_for_individual(name)

        # Create genome for the replaced individual
        genome = self.make_genome(self.generation.newborns[name])

        # Create genome for the replacement individual
        replacement_generation_name = self.fitting_run.get_previous_genetic_generation_name(self.generation.name)
        #print(self.generation.name, replacement_generation_name)
        population = self.fitting_run.get_population_for_generation(replacement_generation_name)
        #replacement_genome = self.make_genome(self.population[replacement_name])
        # OR:
        replacement_genome = self.make_genome(population[replacement_name])

        # Create elitism
        return Elitism(index, genome, replacement_genome, score, replacement_score, name, replacement_name)

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
        for elitism in self.elitisms:

            print("#" + str(elitism.index + 1) + " " + fmt.blue + fmt.underlined + "Elitism replacement:" + fmt.reset)
            print("")

            print("NAMES:")
            print("")

            with fmt.print_in_columns(4) as print_row:

                print_row("Replaced name", ":", elitism.name, "[" + self.generation_name + "]")
                print_row("Replacement name", ":", elitism.replacement_name, "[" + self.previous_generation_name + "]")

            print("")
            print("GENOMES:")
            print("")

            replaced_string = tostr(elitism.replaced.genes, delimiter=" ")
            replacement_string = tostr(elitism.replacement.genes, delimiter=" ")

            with fmt.print_in_columns(3) as print_row:

                print_row("Replaced individual", ":", replaced_string)
                print_row("Replacement individual", ":", replacement_string)

            print("")
            print("PARAMETERS:")
            print("")

            # Get actual parameter values
            replaced_string = self.genome_to_parameters_string(elitism.replaced)
            replacement_string = self.genome_to_parameters_string(elitism.replacement)

            # Print parameters
            with fmt.print_in_columns(3) as print_row:
                print_row("Replaced individual", ":", "[" + replaced_string + "]")
                print_row("Replacement individual", ":", "[" + replacement_string + "]")

            print("")
            print("SCORES:")
            print("")

            # Get the scores
            score = elitism.replaced_score
            replacement_score = elitism.replacement_score

            if replacement_score > score: log.warning("Replacement score is higher than original score")

            # Print scores
            with fmt.print_in_columns(3) as print_row:
                print_row("Replaced score", ":", str(score))
                print_row("Replacement score", ":", str(replacement_score))

            # MAKE SURE THAT THE REPLACED INDIVIDUALS ARE PRESENT IN THE DATABASE (AND THE LOWEST SCORE IN THE STATISTICS)

            print("")
            try:
                individual_key = elitism.replacement_name
                #print("looking for", individual_key)
                individual = get_individual(self.database, self.fitting_run_name, self.generation_index_database, individual_key)
                #print(individual)
                # Get the score
                #score = individual["raw"]
                #fitness = individual["fitness"]
                #print(score, fitness)
                print(fmt.green + "Individual found in the database: elitism executed correctly" + fmt.reset)
            except RuntimeError: print(fmt.red + "Individual not found in the database: elitism has not been executed!" + fmt.reset)

            print("")

            #individuals = get_individuals(self.database, self.fitting_run_name, self.generation_index_database)

            scores = get_individuals_scores(self.database, self.fitting_run_name, self.generation_index_database)
            if not sequences.is_minimum(scores, replacement_score): raise print(fmt.red + "Replacement score " + str(replacement_score) + " is not a minimum of the scores: [" + tostr(scores) + "] of the generation" + fmt.reset)

            scores_previous = get_individuals_scores(self.database, self.fitting_run_name, self.previous_generation_index_database)
            if not sequences.is_minimum(scores_previous, replacement_score): print(fmt.red + "Replacement score " + str(replacement_score) + " is not a minimum of the scores: [" + tostr(scores_previous) + "] of the previous generation (parents)" + fmt.reset)

            print("")

    # -----------------------------------------------------------------

    def check_database(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the database ...")

        # TODO: check whether the individuals are sorted on score

        # Load the elitism table
        if fs.is_file(self.generation.elitism_table_path): elitism_table = self.generation.elitism_table
        else:
            elitism_table = None
            log.warning("No elitism data found for generation '" + self.generation_name + "'")

        # Load the scores table
        scores_table = self.generation.chi_squared_table

        # Load the individuals table
        individuals_table = self.generation.individuals_table

        # Get the scores from the database
        scores_database = get_scores_named_individuals(self.database, self.fitting_run_name, self.generation_index_database)

        # Keep track of the number of elitisms
        nelitisms = 0

        # Keep track of the mismatches
        mismatches = []

        # Loop over the individual names
        for name in scores_database:

            # Check whether name in elitism table as replaced individual, in this case don't check the score
            if elitism_table is not None and name in elitism_table.replaced_names:
                nelitisms += 1
                continue

            # Get the score from the database
            score_database = scores_database[name]

            # Get the simulation name
            if individuals_table.has_individual(name):

                # Get simulation name
                simulation_name = individuals_table.get_simulation_name(name)

                # Get the chi squared value
                score = scores_table.chi_squared_for(simulation_name)

                # Check if equal
                equal = np.isclose(score, score_database, rtol=1.e-4)
                if not equal: mismatches.append((score, score_database))

            # Recurrent
            else: log.debug("Individual '" + name + "' is recurrent, so not present in the individuals table of generation '" + self.generation_name + "'")

        # Report
        if len(mismatches) == 0: log.success(self.generation_name + ": OK")
        else:
            log.error(self.generation_name + ": " + str(len(mismatches)) + " mismatch(es):")
            for score, score_database in mismatches: log.error("   " + str(score) + " , " + str(score_database))
        log.info("Number of elitism replacements: " + str(nelitisms))

    # -----------------------------------------------------------------

    def check_statistics(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the statistics table ...")

        # Determine the generation name
        #generation_name = self.get_generation_name(index)

        # Debugging
        #log.debug("Checking generation '" + generation_name + "' ...")

        # Get the best score from the statistics
        statistics_index = index + 1
        score_statistics = get_best_score_for_generation(self.statistics_path, self.fitting_run_name, statistics_index, minmax="min")

        # Get the best score from the scores table
        score_table = self.get_best_score_for_generation(generation_name)

        # Rel diff
        rel_diff = abs(score_table - score_statistics) / score_statistics

        print(generation_name + ":")
        print("")

        print(" - Best score from statistics: " + tostr(score_statistics))
        print(" - Best score from scores table: " + tostr(score_table))
        print(" - Relative difference: " + tostr(rel_diff) + " (" + tostr(rel_diff * 100) + "%)")

        print("")

    # -----------------------------------------------------------------

    def check_populations(self):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------
