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
import sys
import pexpect

# Import the relevant PTS classes and modules
from pts.core.tools import logging, time
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.core.basics.host import Host
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()
definition.add_required("remote", "string", "the remote host to mount")

# Read the command line arguments
setter = ArgumentConfigurationSetter("mount", "Mount a remote configured in PTS into a local directory")
config = setter.run(definition)

# -----------------------------------------------------------------

# Current working directory
path = fs.cwd()

# Get host
host = Host(config.remote)

#sshfs xxx@nancy.ugent.be: ~/Remotes/nancy -C -o volname=nancy
command = "sshfs " + host.user + "@" + host.name + ": " + path + " -C -o volname=" + host.id

print(command)

# Create the pexpect child instance
child = pexpect.spawn(command, timeout=30)
if host.password is not None:
    child.expect(['password: '])
    child.sendline(host.password)

child.logfile = sys.stdout

# Execute the command and get the output
child.expect(pexpect.EOF, timeout=None)
child.close()

# -----------------------------------------------------------------
