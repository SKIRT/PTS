#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.remote.mounter Contains the RemoteMounter class.

# -----------------------------------------------------------------

# Import standard modules
import sys
import pexpect
import subprocess

# Import the relevant PTS classes and modules
from .host import Host
from .vpn import VPN
from ..tools.logging import log
from ..tools import filesystem as fs
from ..tools import introspection
from .remote import active_keys, add_key

# -----------------------------------------------------------------

# PTS remotes directory
pts_remotes_path = fs.join(introspection.pts_root_dir, "remotes")
if not fs.is_directory(pts_remotes_path): fs.create_directory(pts_remotes_path)

# -----------------------------------------------------------------

class RemoteMounter(object):

    """
    This function ...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(RemoteMounter, self).__init__()

        # The VPN connection
        self.vpn = None

        # The mount paths for the mounted remotes
        self.mount_paths = dict()

        # Check which remotes are mounted
        for path, name in fs.directories_in_path(pts_remotes_path, returns=["path", "name"]):

            try:
                # If empty directory
                if fs.is_empty(path):
                    fs.remove_directory(path)
                    continue
            except OSError: pass

            # Add the path to the dictionary
            self.mount_paths[name] = path

    # -----------------------------------------------------------------

    def is_mounted(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        return host_id in self.mount_paths

    # -----------------------------------------------------------------

    def connect_to_vpn(self, host):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Connecting to vpn service '" + host.vpn.service + "' ...")

        # Connect to the VPN service
        self.vpn = VPN(host.vpn.service)
        self.vpn.connect(host.vpn.user, host.vpn.password, host.vpn.secret, host.vpn.prompt_time_delay)

    # -----------------------------------------------------------------

    def mount(self, host_id):

        """
        This function ...
        :param host_id
        :return:
        """

        # Check if already mounted
        if self.is_mounted(host_id):
            log.warning("The remote host '" + host_id + "' is already mounted")
            return self.mount_paths[host_id]

        # Get host
        host = Host(host_id)

        # If a VPN connection is required for the remote host
        if host.requires_vpn: self.connect_to_vpn(host)

        # Check if key is active
        if host.key is not None:
            if host.key not in active_keys(): add_key(host.key)

        # Create directory for remote
        mount_path = fs.join(pts_remotes_path, host.id)
        #if not fs.is_directory(path): fs.create_directory(path)
        fs.create_directory(mount_path)

        # Debug flags for sshfs
        debug_flags = "-f -d " if log.is_debug() else ""

        # e.g. sshfs xxx@nancy.ugent.be: ~/PTS/remotes/nancy -C -o volname=nancy
        command = "sshfs " + debug_flags + host.user + "@" + host.name + ": " + mount_path + " -C -o volname=" + host.id

        # Create the pexpect child instance
        child = pexpect.spawn(command, timeout=30)
        if host.password is not None:
            child.expect(['password: '])
            child.sendline(host.password)

        child.logfile = sys.stdout

        # Execute the command and get the output
        child.expect(pexpect.EOF, timeout=None)
        child.close()

        # Set the path
        self.mount_paths[host_id] = mount_path

        # Return the path
        return mount_path

    # -----------------------------------------------------------------

    def unmount(self, host_id):

        """
        This function ...
        :return:
        """

        # If not yet mounted
        if not self.is_mounted(host_id):
            log.warning(host_id + " was not mounted")
            return

        # Get host
        host = Host(host_id)

        # Create directory for remote
        path = fs.join(pts_remotes_path, host.id)

        # Unmount
        subprocess.call(["umount", path])

        # Remove the directory for the remote
        fs.remove_directory(path)

    # -----------------------------------------------------------------

    def __del__(self):

        """
        This function ...
        :return:
        """

        # Disconnect from the VPN service if necessary
        #if self.vpn is not None: self.vpn.disconnect()

        # NO, don't do this: we want to keep the connection and the mounted remote open

        pass

# -----------------------------------------------------------------
