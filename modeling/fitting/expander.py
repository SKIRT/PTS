#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.expander Contains the ParameterExpander class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from collections import OrderedDict, defaultdict

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.basics.log import log
from ...core.tools.utils import lazyproperty, memoize_method
from .tables import WeightsTable
from ...core.tools import filesystem as fs
from ...core.tools import tables
from .modelanalyser import FluxDifferencesTable
from ...core.tools import sequences
from ...core.misc.fluxes import ObservedFluxCalculator
from .tables import ChiSquaredTable, BestParametersTable, ModelProbabilitiesTable, ParameterProbabilitiesTable
from ...core.basics.distribution import Distribution
from ...core.tools import formatting as fmt
from ...core.tools.stringify import tostr
from ...core.tools import nr, numbers
from ...core.plot.distribution import plot_distributions, plot_distribution
from .run import FittingRun
from .tables import GenerationsTable
from ...core.filter.filter import parse_filter
from .generation import GenerationInfo
from .weights import WeightsCalculator
from ...core.launch.manager import SimulationManager

# -----------------------------------------------------------------

class ParameterExpander(FittingComponent):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(ParameterExpander, self).__init__(*args, **kwargs)

        # The generation info
        self.info = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 6. Generate the model parameters
        self.generate_models()

        # 7. Set the paths to the input files
        #if self.needs_input: self.set_input()

        # 8. Adjust the ski template
        #self.adjust_ski()

        # 9. Fill the tables for the current generation
        #self.fill_tables()

        # Launch the models
        self.launch()

        # Show
        self.show()

        # Write
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ParameterExpander, self).setup(**kwargs)

        # Load the generation info
        self.info = self.fitting_run.get_generation_info(self.config.generation)

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_run(self):

        """
        This function ...
        :return:
        """

        return self.load_fitting_run(self.config.run)

    # -----------------------------------------------------------------

    @property
    def parameter_ranges(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.free_parameter_ranges

    # -----------------------------------------------------------------

    @property
    def parameter_labels(self):

        """
        This function ...
        :return:
        """

        return self.parameter_ranges.keys()

    # -----------------------------------------------------------------

    @property
    def parameter_units(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.parameter_units

    # -----------------------------------------------------------------

    def get_parameter_unit(self, label):

        """
        Thisf unction ...
        :param label:
        :return:
        """

        return self.parameter_units[label]

    # -----------------------------------------------------------------

    @property
    def initial_parameter_values(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.first_guess_parameter_values

    # -----------------------------------------------------------------

    @lazyproperty
    def generation(self):

        """
        Thisfunction ...
        :return:
        """

        return self.fitting_run.get_generation(self.config.generation)

    # -----------------------------------------------------------------

    @property
    def grid_settings(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.grid_settings

    # -----------------------------------------------------------------

    @lazyproperty
    def parameter_scales(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary for the scales
        scales = dict()

        # Get the scales for each free parameter
        for label in self.parameter_labels:
            key = label + "_scale"
            scales[label] = self.grid_settings[key]

        # Return the scales dict
        return scales

    # -----------------------------------------------------------------

    @property
    def parameters_table(self):

        """
        This function ...
        :return:
        """

        return self.generation.parameters_table

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_parameter_values(self):

        """
        This function ...
        :return:
        """

        return self.parameters_table.unique_parameter_values

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_parameter_values_scalar(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary
        values_scalar = OrderedDict()

        # Loop over the parameters
        for label in self.unique_parameter_values:
            value = self.unique_parameter_values[label].to(self.get_parameter_unit(label)).value
            values_scalar[label] = value

        # Return the scalar values
        return values_scalar

    # -----------------------------------------------------------------

    @property
    def chi_squared_table(self):

        """
        This function ...
        :return:
        """

        return self.generation.chi_squared_table

    # -----------------------------------------------------------------

    def generate_models(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the new model parameters ...")

        for label in self.parameter_labels:

            unique_values = self.unique_parameter_values[label]
            print(label, unique_values)

    # -----------------------------------------------------------------

    def launch(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching ...")

        # Create manager
        #self.manager = SimulationManager()

        # Run the manager
        #self.manager.run(assignment=assignment, timing=self.timing_table, memory=self.memory_table)

        # status=status, info_tables=[parameters, chi_squared], remotes=remotes, simulations=simulations)

        # Set the actual number of simulations for this generation
        #self.generation_info.nsimulations = self.nmodels

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

# -----------------------------------------------------------------
