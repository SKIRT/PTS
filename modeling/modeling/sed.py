#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.modeling.sed Contains the SEDModeler class.
#  Perform radiative transfer modeling using genetic algorithms on a simple SED

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ..fitting.configuration import FittingConfigurer
from ..fitting.initialization.sed import SEDFittingInitializer
from .base import ModelerBase
from ..component.sed import get_ski_template, get_observed_sed, get_sed_plot_path
from ...core.basics.range import IntegerRange, QuantityRange
from ...core.basics.configuration import ConfigurationDefinition, InteractiveConfigurationSetter
from ...core.plot.sed import SEDPlotter

# -----------------------------------------------------------------

class SEDModeler(ModelerBase):

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

        # 3. Do the fitting
        self.fit()

        # 4. Writing
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

        # Plot SED
        if "plot_sed" not in self.history: self.plot_sed()

    # -----------------------------------------------------------------

    def plot_sed(self):

        """
        This function ...
        :return:
        """

        # Add an entry to the history
        self.history.add_entry("plot_sed")

        # Create SED plotter
        plotter = SEDPlotter()

        # Add the observed SED
        sed = get_observed_sed(self.modeling_path)
        plotter.add_sed(sed, "Observations")

        # Run the plotter
        plotter.run(output=get_sed_plot_path(self.modeling_path))

        # Mark the end and save the history file
        self.history.mark_end()
        self.history.save()

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

        # Create a definition
        definition = ConfigurationDefinition()
        default_npackages = max(int(1e4), int(ski.packages()/10))
        definition.add_optional("npackages", "positive_integer", "the number of photon packages per wavelength for the initial generation", default=default_npackages)
        definition.add_flag("selfabsorption", "enable dust self-absorption", default=ski.dustselfabsorption())
        definition.add_flag("transient_heating", "enable transient heating", default=ski.transientheating())

        # Add option for the range of the number of wavelengths
        nwavelengths = ski.nwavelengths()
        min_nwavelengths = max(int(0.1 * nwavelengths), 45)
        max_nwavelengths = max(5 * nwavelengths, 5 * min_nwavelengths)
        default_nwavelengths_range = IntegerRange(min_nwavelengths, max_nwavelengths)
        definition.add_optional("nwavelengths_range", "integer_range", "range for the number of wavelengths to vary over the generations", default=default_nwavelengths_range)
        definition.add_optional("ngrids", "positive_integer", "number of wavelength grids to be generated", default=10)
        definition.add_flag("add_emission_lines", "add additional points to the wavelength grids to sample important dust/gas emission lines", default=False)
        default_wavelength_range = QuantityRange(ski.min_wavelength, ski.max_wavelength)
        definition.add_optional("wavelength_range", "quantity_range", "wavelength range for all wavelength grids", default=default_wavelength_range)

        # Create the setter
        setter = InteractiveConfigurationSetter("Initialization of the ski template", add_cwd=False, add_logging=False)
        config = setter.run(definition, prompt_optional=True)

        # Set fixed settings for the ski model
        initializer.config.npackages = config.npackages
        initializer.config.selfabsorption = config.selfabsorption
        initializer.config.transient_heating = config.transient_heating

        # Set options for the wavelength grids
        initializer.config.wg.npoints_range = config.nwavelengths_range
        initializer.config.wg.ngrids = config.ngrids
        initializer.config.wg.add_emission_lines = config.add_emission_lines
        initializer.config.wg.min_wavelength = config.wavelength_range.min
        initializer.config.wg.max_wavelength = config.wavelength_range.max

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
