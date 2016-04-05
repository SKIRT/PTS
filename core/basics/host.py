#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.host Contains the Host class and a few convenience functions

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os

# Import the relevant PTS classes and modules
from ..tools import configuration, inspection, filesystem
from ..tools.logging import log

# -----------------------------------------------------------------

def find_host_ids():

    """
    This function ...
    :return:
    """

    # Create a list to contain the host ids
    ids = []

    # Search for files that define remote host configurations
    hosts_directory = filesystem.join(inspection.pts_user_dir, "hosts")
    if not filesystem.is_directory(hosts_directory): filesystem.create_directory(hosts_directory, recursive=True)

    # If the hosts directory is empty, place a template host configuration file there and exit with an error
    if len([item for item in os.listdir(hosts_directory) if filesystem.is_file(filesystem.join(hosts_directory, item))]) == 0:

        config = configuration.new()
        config.name = "server.institute.com"
        config.user = "user000"
        config.password = None
        config.output_path = "~/DATA/SKIRT"
        config.scheduler = True
        config.mpi_command = "mpirun"
        config.force_process_binding = False
        config.use_hyperthreading = False
        config.modules = ["examplemodule/2016/version2", "examplemodule2/2016/version0.1.3"]
        config.clusters.default = "cluster_a"
        config.clusters.cluster_a.cores = 16
        config.clusters.cluster_b.cores = 32

        config_file_path = filesystem.join(hosts_directory, "template.cfg")
        config_file = open(config_file_path, 'w')
        config.save(config_file)

        log.error("No remote configuration files were found. Placing a template into PTS/user/hosts. Adjust it"
                  " for the remote hosts you want to use before using 'pts launch' or 'pts status'.")
        exit()

    # Loop over the configuration files in the hosts directory
    for filename in os.listdir(hosts_directory):

        # Skip the template configuration file
        if filename == "template.cfg": continue

        # Determine the full path to the host file
        file_path = filesystem.join(hosts_directory, filename)

        # Ignore directories and hidden files
        if filename.startswith(".") or not os.path.isfile(file_path): continue

        # Get the host id for this line
        host_id = filename.split(".")[0]

        # Add the id to the list
        ids.append(host_id)

    # Return the list of host ids
    return ids

# -----------------------------------------------------------------

def has_simulations(host_id):

    """
    This function ...
    :param host_id:
    :return:
    """

    # Check whether there are simulation files corresponding to this host ID
    host_run_dir = filesystem.join(inspection.skirt_run_dir, host_id)

    # If the host run directory does not exist yet, create it
    if not filesystem.is_directory(host_run_dir): filesystem.create_directory(host_run_dir)

    # If there are no simulation files for this host, skip it
    return len([item for item in os.listdir(host_run_dir) if item.endswith(".sim") and not item.startswith(".")]) > 0

# -----------------------------------------------------------------

class Host(object):

    """
    This class ...
    """

    def __init__(self, host_id, cluster=None):

        """
        The constructor ...
        :param host_id:
        :param cluster:
        :return:
        """

        # -- Attributes --

        # The host ID
        self.id = host_id

        # The name of the cluster
        self.cluster_name = None

        ## Read the host configuration file

        # Determine the path to the configuration file for the specified host and check if it is present
        host_file_path = filesystem.join(inspection.pts_user_dir, "hosts", host_id + ".cfg")
        if not os.path.isfile(host_file_path): raise ValueError("The configuration settings for remote host " + host_id + " could not be found in the PTS/user/hosts directory")

        # Open the host configuration file
        config = configuration.open(host_file_path)

        # Set the host ID and cluster name (if a scheduling system is used)
        if config.scheduler:

            # If no cluster name is given, use the default cluster (defined in the configuration file)
            if cluster is None: self.cluster_name = config.clusters.default

            # If a cluster name is given, use that one
            else: self.cluster_name = cluster

        # Set the entries in the configuration object as attributes of this object
        self.name = config.name
        self.user = config.user
        self.password = config.password
        self.output_path = config.output_path
        self.scheduler = config.scheduler
        self.mpi_command = config.mpi_command
        self.force_process_binding = config.force_process_binding
        self.use_hyperthreading = config.use_hyperthreading
        self.modules = config.modules
        self.clusters = config.clusters
        self.vpn = config.vpn

    # -----------------------------------------------------------------

    @property
    def system_name(self):

        """
        This function ...
        :return:
        """

        if self.cluster_name is not None: return self.id + "-" + self.cluster_name
        else: return self.id

    # -----------------------------------------------------------------

    @property
    def requires_vpn(self):

        """
        This function ...
        :return:
        """

        return self.vpn.service is not None

# -----------------------------------------------------------------
