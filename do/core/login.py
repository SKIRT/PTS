#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.login Login to a remote host configured in PTS.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import subprocess
import tempfile

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.core.remote.host import find_host_ids, Host
from pts.core.tools import filesystem as fs
from pts.core.tools import introspection
from pts.core.tools import terminal
from pts.core.tools import time

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()
definition.add_required("remote", "string", "remote host to mount", choices=find_host_ids())

# Read the command line arguments
setter = ArgumentConfigurationSetter("login", "Login to a remote host configured in PTS")
config = setter.run(definition)

# -----------------------------------------------------------------

# Check platform
if not introspection.is_macos(): raise RuntimeError("This command only works on MacOS")

# -----------------------------------------------------------------

# Get host
host = Host.from_host_id(config.remote)

# If a VPN connection is required for the remote host
#if host.requires_vpn: self.connect_to_vpn(host)

# Check if key is active
#if host.key is not None:
#    if host.key not in active_keys(): add_key(host.key)

#command = "sshfs " + debug_flags + host.user + "@" + host.name + ": " + mount_path + " -C -o volname=" + host.id

# SSH COMMAND: "ssh -p 2935 -X sjversto@nancy.ugent.be"

ssh_path = "/usr/bin/ssh"

#output = terminal.execute("which ssh")
#print(output)

temp_path = introspection.pts_temp_dir
temp_command_path = fs.join(temp_path, time.unique_name("pts_login") + ".command")

#with tempfile.NamedTemporaryFile(suffix='.command') as f:
with open(temp_command_path, 'w') as f:

    f.write('#!/bin/sh\n')
    #f.write("python bb.py\n")

    #subprocess.call(['open', '-W', f.name])

    #ssh_command = "ssh "
    ssh_arguments = ""
    if host.port is not None: ssh_arguments += "-p " + str(host.port) + " "
    ssh_arguments += "-X "
    ssh_arguments += host.user
    ssh_arguments += "@"
    ssh_arguments += host.name

    #splitted = ssh_command.split(" ")

    ssh_options = "'" + ssh_arguments + "'"

    f.write(ssh_path + " " + ssh_arguments)

    #command = ['open', '-W', '-a', 'Terminal.app', ssh_path, "--args", ssh_options]

    #print(command)

    #command = ['open', '-W', '-a', 'Terminal.app']

# make executable
subprocess.call(["chmod", "+x", temp_command_path])

command = ['open', '-W', temp_command_path]
subprocess.call(command)

fs.remove_file(temp_command_path)

# -----------------------------------------------------------------
