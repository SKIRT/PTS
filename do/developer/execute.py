#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.developer.execute Execute a line on a remote host and check the output.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.remote.remote import Remote

# -----------------------------------------------------------------

definition = ConfigurationDefinition()
definition.add_required("remote", "string", "remote host ID")
definition.add_required("line", "string", "line to be executed on the remote host")
config = parse_arguments("sessions", definition)

# -----------------------------------------------------------------

# Initialize the remote
remote = Remote()
if not remote.setup(config.remote): raise RuntimeError("The remote host '" + config.remote + "' is not available at the moment")

# -----------------------------------------------------------------

print("")
#print("-----------------------------------------------------------------")
print("OUTPUT")
print("-----------------------------------------------------------------")
print("")

# Execute the line, show the output
output = remote.execute(config.line)
for line in output: print(line)

print("")
print("-----------------------------------------------------------------")
print("")

# -----------------------------------------------------------------
