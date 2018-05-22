#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.reinitialization Contains the FittingReinitializer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ..component.galaxy import GalaxyModelingComponent
from ...core.basics.log import log
from ..misc.select import select_from_model_suite, select_from_fitting_context, select_from_analysis_context

# -----------------------------------------------------------------

class FittingReinitializer(FittingComponent, GalaxyModelingComponent):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Call the construcor of the base classes
        #super(FittingReinitializer, self).__init__(*args, **kwargs)
        FittingComponent.__init__(self, no_config=True)
        GalaxyModelingComponent.__init__(self, *args, **kwargs)

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Load
        self.load_model()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        Thisn function ...
        :param kwargs:
        :return:
        """

        # Call the setup fucntion of the base class
        #super(FittingReinitializer, self).setup(**kwargs)
        FittingComponent.setup(self, **kwargs)
        GalaxyModelingComponent.setup(self, **kwargs)

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

    @property
    def from_analysis_run(self):

        """
        This function ...
        :return:
        """

        return self.config.origin == "analysis_run"

    # -----------------------------------------------------------------

    def load_model(self):

        """
        This fucntion ...
        :return:
        """

        # Inform the user
        log.info("Getting the model ...")

        # Load from model
        if self.from_model:

            # Get the model
            model_name = self.prompt_model()

            # Update the analysis info
            #self.analysis_run_info.model_name = model_name
            #self.analysis_run_info.parameter_values = self.parameter_values

        # Prompt for a fitting run
        elif self.from_fitting_run:

            # Get the model
            run_id, generation_name, simulation_name, chi_squared = self.prompt_fitting()

            # Update the analysis info
            #self.analysis_run_info.fitting_run = run_id
            #self.analysis_run_info.generation_name = generation_name
            #self.analysis_run_info.simulation_name = simulation_name
            #self.analysis_run_info.chi_squared = chi_squared
            #self.analysis_run_info.parameter_values = self.parameter_values

            # Set the name of the corresponding model of the model suite
            #self.analysis_run_info.model_name = self.definition.name

        # From an analysis run
        elif self.from_analysis_run:

            # Get the model
            self.prompt_analysis()

        # Invalid
        else: raise ValueError("Invalid value for 'origin'")

    # -----------------------------------------------------------------

    def prompt_model(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the model from the model suite ...")

        # Select the model
        model_name, ski, definition, input_paths, parameter_values = select_from_model_suite(self.model_suite, adapt=self.config.adapt, name=self.config.model_name)

        # Set attributes
        self.ski = ski
        self.definition = definition
        self.parameter_values = parameter_values
        self.input_paths = input_paths

        # Set the 'distance' parameter value, since instruments are still not adapted from the default template.
        # Instruments will only be added to the ski file later, so the parameter value obtained from the 'distance'
        # label in 'select_from_model_suite' is still incorrect
        # Other instrument properties should have been fixed (with the 'fix_labels' function)
        self.parameter_values["distance"] = self.galaxy_distance

        # Return the model name
        return model_name

    # -----------------------------------------------------------------

    def prompt_fitting(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the model from the fitting context ...")

        # Select the model
        run_id, generation_name, simulation_name, fitting_run, chi_squared, ski, input_paths, parameter_values = select_from_fitting_context(self.fitting_context)

        # Load the model definition
        definition = fitting_run.model_definition

        # Set attributes
        self.ski = ski
        self.definition = definition
        self.parameter_values = parameter_values
        self.input_paths = input_paths

        # Return the fitting run ID, generation name, simulation name, and chi_squared
        return run_id, generation_name, simulation_name, chi_squared

    # -----------------------------------------------------------------

    def prompt_analysis(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the model from an analysis run ...")

        # Select
        run_name, analysis_run = select_from_analysis_context(self.analysis_context)

        # TODO: further implement ...

# -----------------------------------------------------------------
