#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.configuration Contains the FittingConfigurer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.tools import filesystem as fs
from ...core.tools import introspection
from ...core.simulation.skifile import LabeledSkiFile
from ...core.tools.logging import log
from ...core.tools import parsing
from .generations import GenerationsTable

# -----------------------------------------------------------------

template_ski_path = fs.join(introspection.pts_dat_dir("modeling"), "ski", "labeled_template.ski")
labels_description_path = fs.join(introspection.pts_dat_dir("modeling"), "ski", "labels_description.dat")
labels_types_path = fs.join(introspection.pts_dat_dir("modeling"), "ski", "labels_types.dat")

# -----------------------------------------------------------------

class FittingConfigurer(FittingComponent):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(FittingConfigurer, self).__init__(config)

        # -- Attributes --

        # The ski file template
        self.ski = None

        # The names of the free parameters
        self.parameters = []

        # The descriptions for the free parameters
        self.descriptions = dict()

        # The types for the free parameters
        self.types = dict()

        # The specified ranges for the free parameters
        self.ranges = dict()

        # The names of the filter names
        self.filter_names = []

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the necessary input
        self.load_input()

        # 3. Prompt for the fitting parameters
        self.prompt_parameters()

        # 4. Prompt for the physical parameter ranges
        self.prompt_ranges()

        # 5. Prompt for the fitting filters
        self.prompt_filters()

        # 6. Adjust the labels of the template ski file
        self.adjust_labels()

        # 7. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(FittingConfigurer, self).setup()

    # -----------------------------------------------------------------

    def load_input(self):

        """
        This function ...
        :return:
        """

        # Load the template ski file
        self.load_template()

        # Load the parameter descriptions
        self.load_descriptions()

        # Load the parameter types
        self.load_types()

    # -----------------------------------------------------------------

    def load_template(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the ski file template ...")

        # Load the labeled ski template file
        self.ski = LabeledSkiFile(template_ski_path)

    # -----------------------------------------------------------------

    def load_descriptions(self):

        """
        This function ...
        :return:
        """

        labels, descriptions = np.genfromtxt(labels_description_path, delimiter=" | ", dtype=str, unpack=True)

        # Set the descriptions
        for label, description in zip(labels, descriptions):
            self.descriptions[label] = description

    # -----------------------------------------------------------------

    def load_types(self):

        """
        This function ...
        :return:
        """

        labels, types, default_ranges = np.genfromtxt(labels_types_path, delimiter=" | ", dtype=str, unpack=True)

        # Set the types
        for label, type, default_range in zip(labels, types, default_ranges):
            self.types[label] = (type, default_range)

    # -----------------------------------------------------------------

    def prompt_parameters(self):

        """
        This function ...
        :return:
        """

        # Get all the labels
        labels = self.ski.labels()

        # Get the choices
        indices = get_choices(labels, "free parameters", self.descriptions)

        # Set the chosen free parameters
        for index in indices: self.parameters.append(labels[index])

    # -----------------------------------------------------------------

    def prompt_ranges(self):

        """
        This function ...
        :return:
        """

        # Get the range for each chosen free parameter
        for label in self.parameters:

            # Get the type of the parameter
            parameter_type = self.types[label][0]
            default_range = self.types[label][1]
            if default_range == "--": default_range = None
            range_type = parameter_type + "_range"
            parsing_function = getattr(parsing, range_type)

            log.info("Give a valid physical range for the '" + label + "' parameter:")
            if default_range is not None: log.info("Press ENTER to use the default value (" + str(default_range) + ")")

            # Get the numbers
            answer = raw_input("   : ")

            # Check if default must be used
            if default_range is not None and answer == "": answer = default_range

            # Get the value
            range = parsing_function(answer)

            # Set the range
            self.ranges[label] = range

    # -----------------------------------------------------------------

    def prompt_filters(self):

        """
        This function ...
        :return:
        """

        # Get all the possible filter names
        filter_names = self.observed_filter_names

        # Get the choices
        indices = get_choices(filter_names, "filters")

        # Set the chosen filter names
        for index in indices: self.filter_names.append(filter_names[index])

    # -----------------------------------------------------------------

    def adjust_labels(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adjusting labels of the free parameters in the template ski file ...")

        # Loop over the labels currently in the template ski file
        for label in self.ski.labels():

            # If the label is in the list of free parameter labels, skip (don't remove)
            if label in self.parameters: continue

            # Otherwise, delabel
            self.ski.delabel(label)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the fit parameters
        self.write_parameters()

        # Write the reference filter names
        self.write_filter_names()

        # Write the ski file template
        self.write_ski()

        # Write generations table
        self.write_generations_table()

    # -----------------------------------------------------------------

    def write_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the free parameter labels ...")

        # Write
        with open(self.free_parameters_path, 'w') as f:
            for label in self.parameters: print(label + " | " + str(self.ranges[label].min) + " | " + str(self.ranges[label].max) + " | " + self.descriptions[label], file=f)

    # -----------------------------------------------------------------

    def write_filter_names(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the fitting filter names ...")

        # Write
        with open(self.fitting_filters_path, 'w') as f:
            for name in self.filter_names: print(name, file=f)

    # -----------------------------------------------------------------

    def write_ski(self):

        """
        This function ...
        :return:
        """

        # Save the ski file template
        self.ski.saveto(self.template_ski_path)

    # -----------------------------------------------------------------

    def write_generations_table(self):

        """
        This function ...
        :return:
        """

        # Initialize the generations table
        generations_table = GenerationsTable.initialize(self.parameters)

        # Save the generations table
        generations_table.saveto(self.generations_table_path)

# -----------------------------------------------------------------

def get_choices(options, feature, descriptions=None):

    """
    This function ...
    :param options:
    :param feature:
    :param descriptions:
    :return:
    """

    log.info("Possible " + feature + ":")
    for index, label in enumerate(options):
        description = "  " + descriptions[label] if descriptions is not None else ""
        log.info(" - [" + str(index) + "] " + label + description)
    log.info("")

    log.info("Give the numbers of the " + feature + " that should be used for the fitting:")

    # Get the numbers
    answer = raw_input("   : ")
    indices = parsing.integer_list(answer)

    # Return the chosen indices
    return indices

# -----------------------------------------------------------------
