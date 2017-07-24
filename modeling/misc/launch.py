#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.misc.launch Contains the GalaxyModelingComponent class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ..component.galaxy import GalaxyModelingComponent
from .select import select_from_model_suite, select_from_fitting_context
from ...core.tools.logging import log

# -----------------------------------------------------------------

class ModelLauncher(GalaxyModelingComponent):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ModelLauncher, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup class
        self.setup(**kwargs)

        # 2. Get the model
        self.get_model()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ModelLauncher, self).setup(**kwargs)

    # -----------------------------------------------------------------

    @property
    def from_model(self):

        """
        This function ...
        :return:
        """

        return self.config.origin == "model"

    # -----------------------------------------------------------------

    @property
    def from_fitting_run(self):

        """
        This function ...
        :return:
        """

        return self.config.origin == "fitting_run"

    # -----------------------------------------------------------------

    def get_model(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Getting the model ...")

        # Load from model
        if self.from_model: self.prompt_model()

        # Prompt for a fitting run
        elif self.from_fitting_run: self.prompt_fitting()

        # Invalid
        else: raise ValueError("Invalid value for 'origin'")

        # Show the model parameters
        print("")
        print("Model parameter values:")
        print("")
        for label in self.parameter_values: print(" - " + label + ": " + tostr(self.parameter_values[label]))
        print("")

    # -----------------------------------------------------------------

    def prompt_model(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the model from the model suite ...")

        # NEW
        model_name, ski, definition, input_paths, parameter_values = select_from_model_suite(self.model_suite)

        self.ski = ski
        self.definition = definition
        self.parameter_values = parameter_values
        self.input_paths = input_paths

    # -----------------------------------------------------------------

    def prompt_fitting(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the model from the fitting context ...")

        # NEW
        run_id, generation_name, simulation_name, fitting_run, chi_squared, ski, input_paths, parameter_values = select_from_fitting_context(self.fitting_context)

        # Set
        self.parameter_values = parameter_values
        self.ski = ski
        self.input_paths = input_paths

# -----------------------------------------------------------------
