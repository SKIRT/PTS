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
from ..basics.map import Map
from ..tools import types

# -----------------------------------------------------------------

protocols = ["ssh", "smb"]

# -----------------------------------------------------------------

def find_host_ids(schedulers=None, protocol="ssh"):

    """
    This function ...
    :param schedulers: None, True or False
    :param protocol: ssh, smb or combination
    :return:
    """

    # Set protocol as list
    if types.is_sequence(protocol): pass
    elif types.is_string_type(protocol): protocol = [protocol]
    else: raise ValueError("Invalid option for 'protocol'")

    # Create a list to contain the host ids
    ids = []

    # Search for files that define remote host configurations
    hosts_directory = fs.join(introspection.pts_user_dir, "hosts")
    if not fs.is_directory(hosts_directory): fs.create_directory(hosts_directory, recursive=True)

    # Loop over the configuration files in the hosts directory
    for file_path, host_id in fs.files_in_path(hosts_directory, extension="cfg", returns=["path", "name"]):

        #host = None
        host = load_host(host_id)

        # If schedulers is specified (False or True)
        if schedulers is not None:
            #if host is None: host = load_host(host_id)
            if schedulers != host.scheduler: continue

        # Check whether this host's protocol is part of those requested
        if host.protocol not in protocol: continue

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
    #return len([item for item in os.listdir(host_run_dir) if item.endswith(".sim") and not item.startswith(".")]) > 0
    return fs.nfiles_in_path(host_run_dir, extension="sim") > 0

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
    #return len([item for item in os.listdir(host_run_dir) if item.endswith(".task") and not item.startswith(".")]) > 0
    return len([item for item in os.listdir(host_run_dir) if item.endswith(".task") and not item.startswith(".")]) > 0

# -----------------------------------------------------------------

def load_host(host_id, clustername=None):

    """
    This function ...
    :param host_id:
    :param clustername:
    :return:
    """

    ## Read the host configuration file

    # Determine the path to the configuration file for the specified host and check if it is present
    host_file_path = fs.join(introspection.pts_user_dir, "hosts", host_id + ".cfg")
    if not os.path.isfile(host_file_path): raise ValueError(
        "The configuration settings for remote host " + host_id + " could not be found in the PTS/user/hosts directory")

    # Create and return the host object
    return Host.from_file(host_file_path, clustername=clustername)

# -----------------------------------------------------------------

class Host(object):

    """
    This class ...
    """

    def __init__(self, host_id, **kwargs):

        """
        The constructor ...
        :param host_id:
        :param **kwargs
        :return:
        """

        # -- Attributes --

        # The host ID
        self.id = host_id

        # The name of the cluster
        self.cluster_name = kwargs.pop("cluster_name", None)

        # Set the entries in the configuration object as attributes of this object
        self.name = kwargs.pop("name")
        self.user = kwargs.pop("user")
        self.password = kwargs.pop("password")
        self.scheduler = kwargs.pop("scheduler", False)
        self.scratch_path = kwargs.pop("scratch_path", None)
        self.output_path = kwargs.pop("output_path", None)
        self.mount_point = kwargs.pop("mount_point", None)
        self.mpi_command = kwargs.pop("mpi_command", None)
        self.force_process_binding = kwargs.pop("force_process_binding", False)
        self.use_hyperthreading = kwargs.pop("use_hyperthreading", False)
        self.maximum_walltime = kwargs.pop("maximal_walltime", None)
        self.preferred_walltime = kwargs.pop("preferred_waltime", None)
        self.clusters = kwargs.pop("clusters", Map()) # mapping
        self.vpn = kwargs.pop("vpn", Map()) # mapping
        self.port = kwargs.pop("port", None)
        self.key = kwargs.pop("key", None)
        self.key_password = kwargs.pop("key_password", None)
        self.quota_command = kwargs.pop("quota_command", None)
        self.protocol = kwargs.pop("protocol", "ssh")

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, filepath, clustername=None):

        """
        This function ...
        :param filepath:
        :param clustername:
        :return:
        """

        # Determine the host ID
        host_id = fs.strip_extension(fs.name(filepath))

        # Open the host configuration file
        config = Configuration.from_file(filepath)

        # Set the host ID and cluster name (if a scheduling system is used)
        if config.scheduler:

            # If no cluster name is given, use the default cluster (defined in the configuration file)
            if clustername is None: cluster_name = config.clusters.default

            # If a cluster name is given, use that one
            else: cluster_name = clustername

        # Not a scheduling system
        else: cluster_name = None

        # Set cluster name
        config["cluster_name"] = cluster_name

        # Create and return
        return cls(host_id, **config)

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
