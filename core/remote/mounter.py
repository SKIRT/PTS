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
from .host import Host, load_host
from .vpn import VPN
from ..tools.logging import log
from ..tools import filesystem as fs
from ..tools import introspection
from .remote import active_keys, add_key
from ..tools import types

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

            # Check if this is an existing mount point
            if fs.is_mount_point(path):

                # Try looping over the files in the mounted directory
                try:
                    for _ in fs.files_in_path(path): pass
                    # If this doesn't give an error, this is a valid mount point
                except OSError: # mount point in limbo (connection has been lost but not unmounted)
                    # unmount
                    subprocess.call(["umount", path])
                    # Remove the directory
                    fs.remove_directory(path)
                    # Continue
                    continue

            # Just a directory
            elif fs.is_directory(path):
                if fs.is_empty(path):
                    fs.remove_directory(path)
                    continue
                else: raise RuntimeError("An error occured: directory '" + path + "' not empty but also not mounted")

            else: continue

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

    def mount(self, host_id, path=None):

        """
        This function ...
        :param host_id
        :param path:
        :return:
        """

        # Get host ID
        if isinstance(host_id, Host): the_host_id = host_id.id
        elif types.is_string_type(host_id): the_host_id = host_id
        else: raise ValueError("Invalid value for 'host_id'")

        # Check if already mounted
        if self.is_mounted(the_host_id):
            log.warning("The remote host '" + the_host_id + "' is already mounted")
            return self.mount_paths[the_host_id]

        # Get host
        if isinstance(host_id, Host): host = host_id
        elif types.is_string_type(host_id): host = load_host(host_id)
        else: raise ValueError("Invalid value for 'host_id'")

        # If a VPN connection is required for the remote host
        if host.requires_vpn: self.connect_to_vpn(host)

        # Check if key is active
        if host.key is not None:
            if host.key not in active_keys(): add_key(host.key)

        # Create directory for remote
        if path is not None: mount_path = fs.create_directory_in(path, host.id)
        else: mount_path = fs.create_directory_in(pts_remotes_path, host.id)

        # Debug flags for sshfs
        debug_flags = "-f -d " if log.is_debug() else ""

        # If mount point is defined
        if host.mount_point is not None: mount_point_string = "/" + host.mount_point
        else: mount_point_string = ""

        # e.g. sshfs xxx@nancy.ugent.be: ~/PTS/remotes/nancy -C -o volname=nancy
        # SMB EXAMPLE: mount_smbfs //sjversto:password@files.ugent.be/sjversto/www/users ~/Remotes/WWW
        if host.protocol == "ssh": command = "sshfs " + debug_flags + host.user + "@" + host.name + ":" + mount_point_string + " " + mount_path + " -C -o volname=" + host.id
        elif host.protocol == "smb":
            #fs.remove_directory(mount_path)
            command = "mount_smbfs //" + host.user + "@" + host.name + mount_point_string + " " + mount_path
        else: raise ValueError("Unknown host protocol: " + host.protocol)

        # Debugging
        log.debug("Mounting command: '" + command + "'")

        # Create the pexpect child instance
        child = pexpect.spawn(command, timeout=30)
        if host.password is not None:
            child.expect(['password: ', 'Password for ' + host.name + ": "])
            child.sendline(host.password)

        if log.is_debug(): child.logfile = sys.stdout

        # Execute the command and get the output
        child.expect(pexpect.EOF, timeout=None)
        output = child.before
        child.close()

        lines = output.split("\r\n")
        for line in lines:
            if "error" in line: raise RuntimeError("Something went wrong: " + line)

        if not fs.is_mount_point(mount_path): raise RuntimeError("An error occured during the mounting")

        # Set the path
        self.mount_paths[the_host_id] = mount_path

        # Return the path
        return mount_path

    # -----------------------------------------------------------------

    def unmount(self, host_id):

        """
        This function ...
        :return:
        """

        # Get host ID
        if isinstance(host_id, Host): the_host_id = host_id.id
        elif types.is_string_type(host_id): the_host_id = host_id
        else: raise ValueError("Invalid value for 'host_id'")

        # If not yet mounted
        if not self.is_mounted(the_host_id):
            log.warning(the_host_id + " was not mounted")
            return

        # Get host
        if isinstance(host_id, Host): host = host_id
        elif types.is_string_type(host_id): host = load_host(host_id)
        else: raise ValueError("Invalid value for 'host_id'")

        # Create directory for remote
        path = fs.join(pts_remotes_path, host.id)

        # Unmount
        output = subprocess.check_output(["umount", path])

        # Busy?
        if "busy" in output: raise RuntimeError("Could not unmount: resource is busy")
        else:

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
