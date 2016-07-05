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
from ...core.tools import inspection
from ...core.simulation.skifile import LabeledSkiFile
from ...core.tools.logging import log
from ...core.tools import parsing

# -----------------------------------------------------------------

template_ski_path = fs.join(inspection.pts_dat_dir("modeling"), "ski", "labeled_template.ski")
labels_description_path = fs.join(inspection.pts_dat_dir("modeling"), "ski", "labels_description.dat")

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

        # 4. Prompt for the fitting filters
        self.prompt_filters()

        # 5. Adjust the labels of the template ski file
        self.adjust_labels()

        # 6. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(FittingComponent, self).setup()

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
            for label in self.parameters: print(label + " | " + self.descriptions[label], file=f)

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
    indices = parsing.int_list(answer)

    # Return the chosen indices
    return indices

# -----------------------------------------------------------------
