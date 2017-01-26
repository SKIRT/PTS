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

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
from ..fitting.explorer import ParameterExplorer
from ..fitting.sedfitting import SEDFitter
from ..component.component import load_modeling_history, get_config_file_path, load_modeling_configuration
from ...core.launch.synchronizer import RemoteSynchronizer
from ...core.remote.remote import is_available
from ...core.prep.deploy import Deployer
from ...core.remote.host import Host

# -----------------------------------------------------------------

class Modeler(Configurable):

    """
    This class ...
    """

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

        # Host ids of the available hosts
        self.available_host_ids = set()

        # The modeling history
        self.history = None

    # -----------------------------------------------------------------

    @property
    def allowed_host_ids(self):

        """
        This function ...
        :return:
        """

        if self.config.local: return []
        elif self.config.remote is not None: return [self.config.remote]
        elif self.modeling_config.host_ids is None: return []
        else: return self.modeling_config.host_ids

    # -----------------------------------------------------------------

    @property
    def allowed_fitting_host_ids(self):

        """
        This function ...
        :return:
        """

        if self.config.fitting_local: return []
        elif self.config.fitting_remote is not None: return [self.config.fitting_remote]
        elif self.modeling_config.fitting_host_ids is None: return []
        else: return self.modeling_config.fitting_host_ids

    # -----------------------------------------------------------------

    @property
    def has_allowed_host_ids(self):

        """
        This function ...
        :return:
        """

        return len(self.allowed_host_ids) > 0

    # -----------------------------------------------------------------

    @property
    def has_allowed_fitting_host_ids(self):

        """
        This function ...
        :return:
        """

        return len(self.allowed_fitting_host_ids) > 0

    # -----------------------------------------------------------------

    @property
    def host_id(self):

        """
        This function returns the ID of the host used for the computationaly heavy stuff
        :return:
        """

        if self.allowed_host_ids is None: return None

        # Loop over the preferred hosts
        for host_id in self.allowed_host_ids:
            if not self.config.check_hosts or host_id in self.available_host_ids: return host_id

        # No host avilable
        return None

    # -----------------------------------------------------------------

    @property
    def fitting_host_ids(self):

        """
        This function ...
        :return:
        """

        host_ids = []

        # Loop over the specified hosts
        for host_id in self.allowed_fitting_host_ids:
            if not self.config.check_hosts or host_id in self.available_host_ids: host_ids.append(host_id)

        # Return the list of host IDs
        return host_ids

    # -----------------------------------------------------------------

    @property
    def all_host_ids(self):

        """
        This function ...
        :return:
        """

        host_ids = set()

        # Add main host ID
        if self.host_id is not None: host_ids.add(self.host_id)

        # Add fitting host ids, if they are available
        for host_id in self.fitting_host_ids: host_ids.add(host_id)

        # Return the list of host IDs
        return list(host_ids)

    # -----------------------------------------------------------------

    @property
    def fitting_host_ids_non_schedulers(self):

        """
        This function ...
        :return:
        """

        host_ids = []
        for host_id in self.fitting_host_ids:
            host = Host(host_id)
            if host.scheduler: continue
            else: host_ids.append(host_id)
        return host_ids

    # -----------------------------------------------------------------

    @property
    def fitting_host_ids_schedulers(self):

        """
        This function ...
        :return:
        """

        host_ids = []
        for host_id in self.fitting_host_ids:
            host = Host(host_id)
            if not host.scheduler: continue
            else: host_ids.append(host_id)
        return host_ids

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
        for host_id in self.allowed_host_ids:
            if is_available(host_id):
                log.debug("Host '" + host_id + "' is available")
                self.available_host_ids.add(host_id)
            else: log.debug("Host '" + host_id + "' is not available")

        # Find available hosts from fitting.host_ids list
        for host_id in self.allowed_fitting_host_ids:
            if self.allowed_host_ids is not None and host_id in self.allowed_host_ids: continue
            if is_available(host_id):
                log.debug("Host '" + host_id + "' is available")
                self.available_host_ids.add(host_id)
            else: log.debug("Host '" + host_id + "' is not available")

        # No available host in the list of preferred host ids
        if (self.has_allowed_host_ids or self.has_allowed_fitting_host_ids) and len(self.available_host_ids) == 0:
            raise RuntimeError("None of the preferred hosts are available at this moment")

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
        deployer.config.host_ids = self.all_host_ids

        # Set the host id on which PTS should be installed (on the host for extra computations and the fitting hosts
        # that have a scheduling system to launch the pts run_queue command)
        pts_host_ids = []
        if self.host_id is not None: pts_host_ids.append(self.host_id)
        for host_id in self.fitting_host_ids_schedulers: pts_host_ids.append(host_id)
        deployer.config.pts_on = pts_host_ids

        # Set
        deployer.config.check = self.config.check_versions

        # Run the deployer
        deployer.run()

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
        explorer.config.remotes = self.fitting_host_ids

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
