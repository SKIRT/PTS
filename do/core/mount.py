#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.mount Mount a remote configured in PTS into a local directory.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import platform
import sys
import pexpect
import subprocess

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.core.basics.host import Host
from pts.core.tools import filesystem as fs
from pts.core.tools import introspection
from pts.core.tools.logging import log

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()
definition.add_required("remote", "string", "the remote host to mount")

# Read the command line arguments
setter = ArgumentConfigurationSetter("mount", "Mount a remote configured in PTS into a local directory")
config = setter.run(definition)

# -----------------------------------------------------------------

def open_file(path):
    if platform.system() == "Darwin": subprocess.Popen(["open", path])
    else: subprocess.Popen(["xdg-open", path])

# -----------------------------------------------------------------

# Get host
host = Host(config.remote)

# PTS remotes directory
pts_remotes_path = fs.join(introspection.pts_root_dir, "remotes")
if not fs.is_directory(pts_remotes_path): fs.create_directory(pts_remotes_path)

# Create directory for remote
path = fs.join(pts_remotes_path, host.id)
if not fs.is_directory(path): fs.create_directory(path)

# If not yet mounted
if len(fs.files_in_path(path)) == 0:

    debug_flags = "-f -d " if log.is_debug() else ""

    #e.g. sshfs xxx@nancy.ugent.be: ~/PTS/remotes/nancy -C -o volname=nancy
    command = "sshfs " + debug_flags + host.user + "@" + host.name + ": " + path + " -C -o volname=" + host.id

    # Create the pexpect child instance
    child = pexpect.spawn(command, timeout=30)
    if host.password is not None:
        child.expect(['password: '])
        child.sendline(host.password)

    child.logfile = sys.stdout

    # Execute the command and get the output
    child.expect(pexpect.EOF, timeout=None)
    child.close()

# Open the directory
open_file(path)

# -----------------------------------------------------------------
