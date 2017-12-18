#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.relaunch_simulations_screen Relaunch simulations in a screen from a local script file.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.tools import filesystem as fs
from pts.core.tools import strings
from pts.core.simulation.definition import SingleSimulationDefinition
from pts.core.simulation.arguments import SkirtArguments

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()
definition.add_required("filename", "file_path", "script file path")

# Read the command line arguments
config = parse_arguments("relaunch_simulations_screen", definition, description="Relaunch simulations in a screen from a local script file")

# -----------------------------------------------------------------

# Initiliaze variables
filename = fs.strip_extension(fs.name(config.filename))
#host_id = strings.split_at_last("_")[1]
queue_name = strings.split_at_last(filename, "_")[0]
host_id = None
remote_path = None
screen_name = None

# Get launch info
header = fs.get_header_lines(config.filename)
next_remote_path = False
for line in header:
    if next_remote_path:
        remote_path = strings.unquote(line.strip())
        next_remote_path = False
    elif "running SKIRT on remote host" in line: host_id = strings.unquote(line.split("on remote host ")[1])
    elif "upload this file to the remote filesystem in the following directory" in line:
        next_remote_path = True
    elif line.startswith("screen -S"):
        screen_name = line.split("screen -S ")[1].split()[0]

# -----------------------------------------------------------------

print("queue name", queue_name)
print("host ID", host_id)
print("remote path", remote_path)
print("screen name", screen_name)

# -----------------------------------------------------------------

# Loop over the simulation lines
for line in fs.read_lines(config.filename):
    if line.startswith("#"): continue
    if not line: continue

    # Get SKIRT arguments
    arguments = SkirtArguments.from_command(line)

    print(arguments)

    # TODO: etc...

# -----------------------------------------------------------------
