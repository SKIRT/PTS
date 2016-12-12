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
from ..tools import introspection, network
from ..tools import filesystem as fs
from ..basics.configuration import Configuration

# -----------------------------------------------------------------

class HostDownException(Exception):

    """
    This exception should be raised when connection to a host is not possible because it is temporarily down
    """

# -----------------------------------------------------------------

def find_host_ids(schedulers=None):

    """
    This function ...
    :param schedulers: None, True or False
    :return:
    """

    # Create a list to contain the host ids
    ids = []

    # Search for files that define remote host configurations
    hosts_directory = fs.join(introspection.pts_user_dir, "hosts")
    if not fs.is_directory(hosts_directory): fs.create_directory(hosts_directory, recursive=True)

    # Loop over the configuration files in the hosts directory
    for filename in os.listdir(hosts_directory):

        # Determine the full path to the host file
        file_path = fs.join(hosts_directory, filename)

        # Ignore directories and hidden files
        if filename.startswith(".") or not os.path.isfile(file_path): continue

        # Get the host id for this line
        host_id = filename.split(".")[0]

        # If schedulers is specified (False or True)
        if schedulers is not None:
            host = Host(host_id)
            if schedulers != host.scheduler: continue

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

    # Check whether the SKIRT run directory can be found
    #if introspection.skirt_run_dir is None: raise RuntimeError("The SKIRT run directory could not be located. Missing SKIRT installation?")
    if introspection.skirt_run_dir is None: return False

    # Check whether there are simulation files corresponding to this host ID
    host_run_dir = fs.join(introspection.skirt_run_dir, host_id)

    # If the host run directory does not exist yet, create it
    if not fs.is_directory(host_run_dir): fs.create_directory(host_run_dir)

    # Check whether simulation files are present
    return len([item for item in os.listdir(host_run_dir) if item.endswith(".sim") and not item.startswith(".")]) > 0

# -----------------------------------------------------------------

def has_tasks(host_id):

    """
    This function ...
    :param host_id:
    :return:
    """

    if not fs.is_directory(introspection.pts_run_dir): fs.create_directory(introspection.pts_run_dir)

    # Check whether there are task files corresponding to this host ID
    host_run_dir = fs.join(introspection.pts_run_dir, host_id)

    # If the host run directory does not exist yet, create it
    if not fs.is_directory(host_run_dir): fs.create_directory(host_run_dir)

    # Check whether task files are present
    return len([item for item in os.listdir(host_run_dir) if item.endswith(".task") and not item.startswith(".")]) > 0

# -----------------------------------------------------------------

class Host(object):

    """
    This class ...
    """

    def __init__(self, host_id, clustername=None):

        """
        The constructor ...
        :param host_id:
        :param clustername:
        :return:
        """

        # -- Attributes --

        # The host ID
        self.id = host_id

        # The name of the cluster
        self.cluster_name = None

        ## Read the host configuration file

        # Determine the path to the configuration file for the specified host and check if it is present
        host_file_path = fs.join(introspection.pts_user_dir, "hosts", host_id + ".cfg")
        if not os.path.isfile(host_file_path): raise ValueError("The configuration settings for remote host " + host_id + " could not be found in the PTS/user/hosts directory")

        # Open the host configuration file
        config = Configuration.from_file(host_file_path)

        # Set the host ID and cluster name (if a scheduling system is used)
        if config.scheduler:

            # If no cluster name is given, use the default cluster (defined in the configuration file)
            if clustername is None: self.cluster_name = config.clusters.default

            # If a cluster name is given, use that one
            else: self.cluster_name = clustername

        # Set the entries in the configuration object as attributes of this object
        self.name = config.name
        self.user = config.user
        self.password = config.password
        self.scheduler = config.scheduler
        self.output_path = config.output_path
        self.mpi_command = config.mpi_command
        self.force_process_binding = config.force_process_binding
        self.use_hyperthreading = config.use_hyperthreading
        self.maximum_walltime = config.maximal_walltime
        self.preferred_walltime = config.preferred_waltime
        self.clusters = config.clusters # mapping
        self.vpn = config.vpn # mapping
        self.port = config.port
        self.key = config.key

    # -----------------------------------------------------------------

    @property
    def cluster(self):

        """
        This function ...
        :return:
        """

        return self.clusters[self.cluster_name]

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

        # If a VPN service is specified in the host configuration file
        if self.vpn.service is not None:

            # If the VPN connection is not required when connected to a certain DNS server, check which servers we
            # are connected to
            if self.vpn.not_for_dns_domain is not None:

                # Get the DNS search domains currently connected to
                domains = network.dns_search_domains()

                # Check whether the DNS domain in the host configuration file is in the list of current DNS domains
                # If this is the case, no VPN connection is required. Else, the VPN connection is required.
                if self.vpn.not_for_dns_domain in domains: return False
                else: return True

            # If the VPN connection is always required, regardless of which DNS server we are connected to, return True.
            else: return True

        # If no VPN service is specified in the host configuration file, assume VPN is not required
        else: return False

# -----------------------------------------------------------------
