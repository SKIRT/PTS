#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.modeling.sed Contains the SEDModeler class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ..fitting.configuration import FittingConfigurer
from ..fitting.initialization.sed import SEDFittingInitializer
from ..fitting.component import get_generations_table
from .modeler import Modeler
from ..component.sed import get_ski_template, get_observed_sed
from ...core.basics.range import IntegerRange

# -----------------------------------------------------------------

class SEDModeler(Modeler):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        """

        # Call the constructor of the base class
        super(SEDModeler, self).__init__(config)

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the data
        self.load_data()

        # 8. Do the fitting
        self.fit()

        # 9. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(SEDModeler, self).setup()

    # -----------------------------------------------------------------

    def load_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the input data ...")

        # Get the galaxy properties
        if "fetch_properties" not in self.history: pass

        # Get the SEDs
        if "fetch_seds" not in self.history: pass

    # -----------------------------------------------------------------

    def fit(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fitting radiative transfer models to the data ...")

        # Configure the fitting
        if "configure_fit" not in self.history: self.configure_fit()

        # Initialize the fitting
        if "initialize_fit_sed" not in self.history: self.initialize_fit()

        # Load the generations table
        generations = get_generations_table(self.modeling_path)

        # If some generations have not finished, check the status of and retrieve simulations
        if generations.has_unfinished: self.synchronize()

        # If some generations have finished, fit the SED
        if generations.has_finished: self.fit_sed()

        # IF all generations have finished, explore new generation of models
        if generations.all_finished: self.explore()

    # -----------------------------------------------------------------

    def configure_fit(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Configuring the fitting ...")

        # Create configuration for the FittingConfigurer
        config = dict()

        # Load the ski template, get the free parameters
        ski = get_ski_template(self.config.path)
        free_parameter_names = ski.labels

        # Load the SED, get the fitting filters
        sed = get_observed_sed(self.config.path)
        fitting_filter_names = sed.filter_names()

        # Set free parameters and fitting filters
        config["parameters"] = free_parameter_names
        config["filters"] = fitting_filter_names

        # Create the fitting configurer
        configurer = FittingConfigurer(config)

        # Add an entry to the history
        self.history.add_entry(FittingConfigurer.command_name())

        # Set the working directory
        configurer.config.path = self.modeling_path

        # Run the fitting configurer
        configurer.run()

        # Mark the end and save the history file
        self.history.mark_end()
        self.history.save()

    # -----------------------------------------------------------------

    def initialize_fit(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Initializing the fitting ...")

        # Create the fitting initializer
        initializer = SEDFittingInitializer()

        # Add an entry to the history
        self.history.add_entry(SEDFittingInitializer.command_name())

        # Set the working directory
        initializer.config.path = self.modeling_path

        # Load the current ski template
        ski = get_ski_template(self.modeling_path)

        # Set fixed settings for the ski model
        initializer.config.npackages = ski.packages()
        initializer.config.selfabsorption = ski.dustselfabsorption()
        initializer.config.transient_heating = ski.transientheating()

        # Set options for the wavelength grids
        nwavelengths = ski.nwavelengths()
        min_nwavelengths = max(int(0.1 * nwavelengths), 20)
        max_nwavelengths = nwavelengths
        initializer.config.wg.npoints_range = IntegerRange(min_nwavelengths, max_nwavelengths)
        initializer.config.wg.ngrids = 5
        initializer.config.wg.add_emission_lines = False
        initializer.config.wg.min_wavelength = ski.minwavelength()
        initializer.config.wg.max_wavelength = ski.maxwavelength()

        # OPTIONS FOR THE DUST GRID NOT RELEVANT FOR SED MODELING (YET)

        # Run the fitting initializer
        initializer.run()

        # Mark the end and save the history file
        self.history.mark_end()
        self.history.save()

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

# -----------------------------------------------------------------
