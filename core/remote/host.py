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
from ..basics.composite import SimplePropertyComposite

# -----------------------------------------------------------------

protocols = ["ssh", "smb"]

# -----------------------------------------------------------------

def all_host_ids():

    """
    This function ...
    :param schedulers:
    :return:
    """

    return find_host_ids(protocol=protocols)

# -----------------------------------------------------------------

def ssh_host_ids():

    """
    This function ...
    :return:
    """

    return find_host_ids(protocol="ssh")

# -----------------------------------------------------------------

def smb_host_ids():

    """
    This function ...
    :return:
    """

    return find_host_ids(protocol="smb")

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

        # Load the host
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

def find_hosts(schedulers=None, protocol="ssh", as_dict=False):

    """
    This function ...
    :param schedulers:
    :param protocol:
    :param as_dict:
    :return:
    """

    if as_dict: return {host_id: load_host(host_id) for host_id in find_host_ids(schedulers=schedulers, protocol=protocol)}
    else: return [load_host(host_id) for host_id in find_host_ids(schedulers=schedulers, protocol=protocol)]

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

    # Check if host ID and clustername are defined as one string
    if ":" in host_id:
        if clustername is not None: raise ValueError("Clustername is defined in host name and separately")
        host_id, clustername = host_id.split(":")

    ## Read the host configuration file

    # Determine the path to the configuration file for the specified host and check if it is present
    host_file_path = fs.join(introspection.pts_user_dir, "hosts", host_id + ".cfg")
    if not os.path.isfile(host_file_path): raise ValueError(
        "The configuration settings for remote host " + host_id + " could not be found in the PTS/user/hosts directory")

    # Create and return the host object
    return Host.from_file(host_file_path, clustername=clustername)

# -----------------------------------------------------------------

class Host(SimplePropertyComposite):

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

        # Call the constructor of the base class
        super(Host, self).__init__()

        # -- Attributes --

        # Define properties
        self.add_string_property("id", "host ID")
        self.add_string_property("cluster_name", "name of the cluster to use")
        self.add_string_property("name", "host name (adress)")
        self.add_string_property("user", "username")
        self.add_string_property("password", "user password")
        self.add_boolean_property("scheduler", "uses scheduling system")
        self.add_string_property("scratch_path", "scratch path")
        self.add_string_property("output_path", "output path")
        self.add_string_property("mount_point", "mount point")
        self.add_string_property("mpi_command", "mpi command")
        self.add_boolean_property("force_process_binding", "force process binding to cores")
        self.add_boolean_property("use_hyperthreading", "enable hyperthreading")
        self.add_real_property("maximum_walltime", "maximal allowed walltime")
        self.add_real_property("preferred_walltime", "maximum walltime preferred by user")

        # Clusters: dynamic section: the number of properties in this mapping can vary from remote to remote
        self.add_section("clusters", "clusters", dynamic=True)

        # VPN
        self.add_section("vpn", "vpn settings")
        self.vpn.add_string_property("service", "VPN service name")
        self.vpn.add_string_property("secret", "secret")
        self.vpn.add_real_property("prompt_time_delay", "time delay between showing the prompt (for username and password) and checking for connection")
        self.vpn.add_string_property("not_for_dns_domain", "don't connect to VPN when on this DNS domain")
        self.vpn.add_string_property("user", "username")
        self.vpn.add_string_property("password", "user password")

        self.add_integer_property("port", "port")
        self.add_string_property("key", "key")
        self.add_string_property("key_password", "key password")
        self.add_string_property("quota_command", "user quota command")
        self.add_string_property("protocol", "protocol (ssh or smb)")

        ## SET

        #print("HERE", "clusters" in self.__dict__)

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

    def saveto(self, path, update_path=True):

        """
        This function ...
        :param path:
        :param update_path:
        :return:
        """

        # Create configuration object
        config = Configuration()
        config.name = self.name
        config.user = self.user
        config.password = self.password
        config.scheduler = self.scheduler
        config.scratch_path = self.scratch_path
        config.output_path = self.output_path
        config.mount_point = self.mount_point
        config.mpi_command = self.mpi_command
        config.force_process_binding = self.force_process_binding
        config.use_hyperthreading = self.use_hyperthreading
        config.maximum_walltime = self.maximum_walltime
        config.preferred_walltime = self.preferred_walltime
        config.clusters = self.clusters
        config.vpn = self.vpn
        config.port = self.port
        config.key = self.key
        config.key_password = self.key_password
        config.quota_command = self.quota_command
        config.protocol = self.protocol

        # Write the host configuration
        config.saveto(path)

        # Change path
        if update_path: self._path = path

    # -----------------------------------------------------------------

    @property
    def has_cluster_name(self):

        """
        This function ...
        :return:
        """

        return self.cluster_name is not None

    # -----------------------------------------------------------------

    @property
    def has_cluster(self):

        """
        This function ...
        :return:
        """

        return self.has_cluster_name and self.has_clusters

    # -----------------------------------------------------------------

    @property
    def cluster_names(self):

        """
        This function ...
        :return:
        """

        names = self.clusters.keys()
        names.remove("default")
        return names

    # -----------------------------------------------------------------

    @property
    def nclusters(self):

        """
        This function ...
        :return:
        """

        return len(self.cluster_names)

    # -----------------------------------------------------------------

    @property
    def has_clusters(self):

        """
        This function ...
        :return:
        """

        return self.nclusters > 0

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

    def as_tuple(self):

        """
        This function ...
        :return:
        """

        return (self.id, self.cluster_name)

    # -----------------------------------------------------------------

    def as_string(self):

        """
        This function ...
        :return:
        """

        if self.cluster_name is not None: return self.id + ":" + self.cluster_name
        else: return self.id

    # -----------------------------------------------------------------

    def __hash__(self):

        """
        This function ...
        :return:
        """

        return hash(self.as_tuple())

    # -----------------------------------------------------------------

    def __eq__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        # Compare to string (host ID)
        if types.is_string_type(other):

            if ":" in other: host_id, clustername = other.split(":")
            else: host_id, clustername = other, None

            if clustername is None: return self.id == host_id
            else: return self.id == host_id and self.cluster_name == clustername

        # Compare to other host object
        elif isinstance(other, Host): return self.id == other.id and self.cluster_name == other.cluster_name

        # Invalid input
        else: return False #raise ValueError("Invalid input: must be string or Host") # NO, JUST RETURN FALSE

# -----------------------------------------------------------------
