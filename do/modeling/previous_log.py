#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.previous_log Show previous log file.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import defaultdict

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import verify_modeling_cwd
from pts.modeling.component.component import get_log_file_paths
from pts.core.tools import filesystem as fs
from pts.core.basics.log import log
from pts.core.tools import time
from pts.core.tools import formatting as fmt

# -----------------------------------------------------------------

# Get the modeling commands
modeling_path = verify_modeling_cwd()
filepaths = get_log_file_paths(modeling_path)

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add settings
definition.add_required("match", "string", "(partial) command name")

# Get configuration
config = parse_arguments("previous_config", definition)

# -----------------------------------------------------------------

matches = defaultdict(list)

# -----------------------------------------------------------------

# Loop over the filepaths
for filepath in filepaths:

    # Get the modeling command
    filename = fs.name(filepath)
    command = filename.split("__")[0]

    # Skip other commands
    if not command.startswith(config.match): continue

    # Add filepath
    matches[command].append(filepath)

# -----------------------------------------------------------------

if len(matches) > 1: raise ValueError("Ambigious command: matches are " + ", ".join(matches.keys()))

# -----------------------------------------------------------------

single_command = matches.keys()[0]
log.debug("The matching command is '" + single_command + "'")
if len(matches[single_command]) == 0: raise RuntimeError("Something went wrong")

# -----------------------------------------------------------------

# Get latest
latest_time = None
latest_filepath = None

# Loop over the files
for filepath in matches[single_command]:

    name = fs.strip_extension(fs.name(filepath))
    command_name, datetime = time.get_name_and_time_from_unique_name(name)

    #print(name, command_name, datetime)

    if latest_time is None or latest_time < datetime:
        latest_time = datetime
        latest_filepath = filepath

# -----------------------------------------------------------------

log.info("Showing log from '" + single_command + "' command at " + str(latest_time))

# -----------------------------------------------------------------

# Loop over the lines
for line in fs.read_lines(latest_filepath):

    identifier = line[23:26].strip()

    # Determine type of message
    if identifier == "D": print(fmt.blue + line + fmt.reset)
    elif identifier == "-": print(fmt.green + line + fmt.reset)
    elif identifier == "!": print(fmt.magenta + line + fmt.reset)
    elif identifier == "*": print(fmt.red + line + fmt.reset)
    else: print(line)

# -----------------------------------------------------------------
