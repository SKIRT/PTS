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
from ...core.basics.configurable import Configurable
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
from ..fitting.configuration import FittingConfigurer
from ..fitting.initialization import FittingInitializer
from ..fitting.explorer import ParameterExplorer
from ..fitting.sedfitting import SEDFitter
from ..core.component import load_modeling_history, get_config_file_path, load_modeling_configuration
from ..fitting.component import get_generations_table
from ...core.launch.synchronizer import RemoteSynchronizer
from ...core.remote.remote import is_available
from ...core.prep.deploy import Deployer

# -----------------------------------------------------------------

# Set the filters for which the data can be used for fitting (data that can be trusted well enough)
fitting_filter_names = ["GALEX FUV", "GALEX NUV", "SDSS u", "SDSS g", "SDSS r", "SDSS i", "SDSS z", "WISE W1",
                        "IRAC I1", "IRAC I2", "WISE W2", "IRAC I3", "IRAC I4", "WISE W3", "WISE W4", "MIPS 24mu",
                        "Pacs blue", "Pacs red", "SPIRE PSW", "SPIRE PMW", "SPIRE PLW"]

# -----------------------------------------------------------------

class SEDModeler(Configurable):

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

        # The path to the modeling directory
        self.modeling_path = None

        # The modeling configuration
        self.modeling_config = None

        # Host ids of the available hosts
        self.available_host_ids = set()

        # The modeling history
        self.history = None

    # -----------------------------------------------------------------

    @property
    def used_host_ids(self):

        """
        This function ...
        :return:
        """

        host_ids = set()

        # Add main host ID
        host_ids.add(self.host_id)

        # Add fitting host ids, if they are available
        for host_id in self.modeling_config.fitting_host_ids:
            if host_id in self.available_host_ids: host_ids.add(host_id)

        # Return the list of host IDs
        return list(host_ids)

    # -----------------------------------------------------------------

    @property
    def host_id(self):

        """
        This function ...
        :return:
        """

        # Loop over the preferred hosts
        for host_id in self.modeling_config.host_ids:
            if host_id in self.available_host_ids: return host_id

        # No host avilable
        return None

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

        # Set the path to the modeling directory
        self.modeling_path = self.config.path

        # Check for the presence of the configuration file
        if not fs.is_file(get_config_file_path(self.modeling_path)): raise ValueError("The current working directory is not a radiative transfer modeling directory (the configuration file is missing)")
        else: self.modeling_config = load_modeling_configuration(self.modeling_path)

        # Find which hosts are available
        if self.config.check_hosts: self.find_available_hosts()

        # Load the modeling history
        self.history = load_modeling_history(self.modeling_path)

        # Deploy SKIRT and PTS
        if self.config.deploy: self.deploy()

    # -----------------------------------------------------------------

    def find_available_hosts(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Finding available hosts ...")

        # Find available hosts from host_ids list
        for host_id in self.modeling_config.host_ids:
            if is_available(host_id):
                log.debug("Host '" + host_id + "' is available")
                self.available_host_ids.add(host_id)
            else: log.debug("Host '" + host_id + "' is not available")

        # Find available hosts from fitting.host_ids list
        for host_id in self.modeling_config.fitting_host_ids:
            if host_id in self.modeling_config.host_ids: continue
            if is_available(host_id):
                log.debug("Host '" + host_id + "' is available")
                self.available_host_ids.add(host_id)
            else: log.debug("Host '" + host_id + "' is not available")

        # No available host in the list of preferred host ids
        if len(self.available_host_ids) == 0: raise RuntimeError("None of the preferred hosts are available at this moment")

    # -----------------------------------------------------------------

    def deploy(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Deploying SKIRT and PTS ...")

        # Create the deployer
        deployer = Deployer()

        # Set the host ids
        deployer.config.host_ids = self.used_host_ids

        # Set the host id on which PTS should be installed
        deployer.config.pts_on = [self.host_id]

        # Set
        deployer.config.check = self.config.check_versions

        # Run the deployer
        deployer.run()

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
        if "initialize_fit" not in self.history: self.initialize_fit()

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

        # Create configuration
        config = dict()

        # Set free parameters
        config["parameters"] = free_parameters[self.modeling_config.method]
        config["ranges"] = free_parameter_ranges[self.modeling_config.method]
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
        initializer = FittingInitializer()

        # Add an entry to the history
        self.history.add_entry(FittingInitializer.command_name())

        # Set the working directory
        initializer.config.path = self.modeling_path

        # Run the fitting initializer
        initializer.run()

        # Mark the end and save the history file
        self.history.mark_end()
        self.history.save()

    # -----------------------------------------------------------------

    def synchronize(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Synchronizing with the remotes (retrieving and analysing finished models) ...")

        # Create the remote synchronizer
        synchronizer = RemoteSynchronizer()

        # Set the host IDs
        synchronizer.config.host_ids = self.modeling_config.fitting_host_ids

        # Run the remote synchronizer
        synchronizer.run()

    # -----------------------------------------------------------------

    def fit_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fitting the SED to the finished generations ...")

        # Create the SED fitter
        fitter = SEDFitter()

        # Add an entry to the history
        self.history.add_entry(SEDFitter.command_name())

        # Run the fitter
        fitter.run()

        # Mark the end and save the history file
        self.history.mark_end()

    # -----------------------------------------------------------------

    def explore(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Exploring the parameter space ...")

        # Create the parameter explorer
        explorer = ParameterExplorer()

        # Add an entry to the history
        self.history.add_entry(ParameterExplorer.command_name())

        # Set the working directory
        explorer.config.path = self.modeling_path

        # Run the parameter explorer
        explorer.run()

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
