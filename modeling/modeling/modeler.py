#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.modeler Contains the GalaxyModeler class, which runs the radiative transfer modelling procedure
#  for a certain galaxy.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta, abstractmethod

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
from ..fitting.explorer import ParameterExplorer
from ..fitting.sedfitting import SEDFitter
from ..component.component import load_modeling_history, get_config_file_path, load_modeling_configuration
from ...core.launch.synchronizer import RemoteSynchronizer
from ...core.prep.deploy import Deployer
from ..fitting.component import get_generations_table, get_ngenerations, has_unevaluated_generations, has_unfinished_generations
from ...core.remote.moderator import PlatformModerator

# -----------------------------------------------------------------

def repeat(target, ntimes):

    """
    This function ...
    :param target:
    :param ntimes:
    :return:
    """

    for _ in range(ntimes): target()

# -----------------------------------------------------------------

class Modeler(Configurable):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        """

        # Call the constructor of the base class
        super(Modeler, self).__init__(config)

        # The path to the modeling directory
        self.modeling_path = None

        # The modeling configuration
        self.modeling_config = None

        # Platform moderator
        self.moderator = None

        # The modeling history
        self.history = None

    # -----------------------------------------------------------------

    @property
    def configured_fitting_host_ids(self):

        """
        This function ...
        :return:
        """

        if self.modeling_config.fitting_host_ids is None: return []
        else: return self.modeling_config.fitting_host_ids

    # -----------------------------------------------------------------

    @property
    def has_configured_fitting_host_ids(self):

        """
        This function ...
        :return:
        """

        return len(self.configured_fitting_host_ids) > 0

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(Modeler, self).setup(**kwargs)

        # Set the path to the modeling directory
        self.modeling_path = self.config.path

        # Check for the presence of the configuration file
        if not fs.is_file(get_config_file_path(self.modeling_path)): raise ValueError("The current working directory is not a radiative transfer modeling directory (the configuration file is missing)")
        else: self.modeling_config = load_modeling_configuration(self.modeling_path)

        # Set execution platforms
        self.set_platforms()

        # Check the number of generations
        if self.config.ngenerations > 1 and self.moderator.any_remote: raise ValueError("When remote execution is enabled, the number of generations per run can only be one")

        # Load the modeling history
        self.history = load_modeling_history(self.modeling_path)

        # Deploy SKIRT and PTS
        if self.config.deploy: self.deploy()

    # -----------------------------------------------------------------

    def set_platforms(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Determining execution platforms ...")

        # Setup the platform moderator
        self.moderator = PlatformModerator()

        # Set platform(s) for fitting (simulations)
        if self.config.fitting_local: self.moderator.add_local("fitting")
        elif self.config.fitting_remotes is not None: self.moderator.add_ensemble("fitting", self.config.fitting_remotes)
        elif self.modeling_config.fitting_host_ids is None: self.moderator.add_local("fitting")
        else: self.moderator.add_ensemble("fitting", self.modeling_config.fitting_host_ids)

        # Other computations
        if self.config.local: self.moderator.add_local("other")
        elif self.config.remotes is not None: self.moderator.add_single("other", self.config.remotes)
        elif self.modeling_config.host_ids is None: self.moderator.add_local("other")
        else: self.moderator.add_single("other", self.modeling_config.host_ids)

        # Run the moderator
        self.moderator.run()

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
        deployer.config.host_ids = self.moderator.all_host_ids

        # Set the host id on which PTS should be installed (on the host for extra computations and the fitting hosts
        # that have a scheduling system to launch the pts run_queue command)
        deployer.config.pts_on = self.moderator.all_host_ids

        # Set
        deployer.config.check = self.config.check_versions

        # Run the deployer
        deployer.run()

    # -----------------------------------------------------------------

    def fit(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fitting radiative transfer models ...")

        # Configure the fitting
        if not self.history.has_configured_fit: self.configure_fit()

        # Initialize the fitting
        if not self.history.has_initialized_fit: self.initialize_fit()

        # If finishing the generation is requested
        if self.config.finish: self.finish()

        # Advance the fitting with a new generations
        else: repeat(self.advance, self.config.ngenerations)

    # -----------------------------------------------------------------

    @abstractmethod
    def configure_fit(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractmethod
    def initialize_fit(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def advance(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Advancing the fitting with a new generation ...")

        # Load the generations table
        generations = get_generations_table(self.modeling_path)

        # Debugging
        if generations.last_generation_name is not None: log.debug("Previous generation: " + generations.last_generation_name)
        else: log.debug("Will be launching the initial generation ...")

        # If some generations have not finished, check the status of and retrieve simulations
        if generations.has_unfinished and self.has_configured_fitting_host_ids: self.synchronize()

        # If some generations have finished, fit the SED
        if generations.has_finished and has_unevaluated_generations(self.modeling_path): self.fit_sed()

        # If all generations have finished, explore new generation of models
        if generations.all_finished: self.explore()

        # Do SED fitting after the exploration step if it has been performed locally (simulations are done, evaluation can be done directly)
        if self.moderator.single_is_local("fitting"): self.finish()

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

        # Set the remote host IDs
        explorer.config.remotes = self.moderator.host_ids_for_ensemble("fitting")

        # Set the number of simulations per generation
        if self.config.nsimulations is not None: explorer.config.nsimulations = self.config.nsimulations

        # Set other settings
        explorer.config.npackages_factor = self.config.npackages_factor
        explorer.config.increase_npackages = self.config.increase_npackages
        explorer.config.refine_wavelengths = self.config.refine_wavelengths
        explorer.config.refine_dust = self.config.refine_dust
        explorer.config.selfabsorption = self.config.selfabsorption
        explorer.config.transient_heating = self.config.transient_heating

        # Run the parameter explorer
        explorer.run()

        # Mark the end and save the history file
        self.history.mark_end()
        self.history.save()

    # -----------------------------------------------------------------

    def finish(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Evaluating last generation ...")

        # Check the current number of generations
        current_ngenerations = get_ngenerations(self.modeling_path)
        if current_ngenerations <= 1: raise RuntimeError("Need at least one generation after the initial generation to finish the fitting")

        # Check if there are unfinished generations
        has_unfinished = has_unfinished_generations(self.modeling_path)
        if has_unfinished: log.warning("There are unfinished generations, but evaluting finished simulations anyway ...")

        # Check if there are unevaluated generations
        if not has_unevaluated_generations(self.modeling_path):
            log.success("All generations have already been evaluated")

        # Do the SED fitting step
        self.fit_sed()

        # Success
        if not has_unfinished: log.success("Succesfully evaluted all generations")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

# -----------------------------------------------------------------
