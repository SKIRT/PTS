#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.show_simulation_log Show the log output of a remote SKIRT simulation.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.core.basics.host import find_host_ids
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs
from pts.core.basics.remote import Remote
from pts.core.tools import introspection
from pts.core.simulation.simulation import RemoteSimulation

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add required
definition.add_required("remote", "string", "name of the remote host", choices=find_host_ids())
definition.add_required("id", "positive_integer", "simulation ID")

# -----------------------------------------------------------------

# Parse the arguments into a configuration
setter = ArgumentConfigurationSetter("show_simulation_log", "Show the log output of a remote SKIRT simulation")
config = setter.run(definition)

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = fs.join(fs.cwd(), time.unique_name("show_simulation_log") + ".txt") if config.report else None

# Determine the log level
level = "DEBUG" if config.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting show_simulation_log ...")

# -----------------------------------------------------------------

# Create and setup the remote
remote = Remote()
remote.setup(config.remote)

# Open the simulation
simulation_path = fs.join(introspection.skirt_run_dir, config.remote, str(config.id) + ".sim")
simulation = RemoteSimulation.from_file(simulation_path)

# The name of the ski file (the simulation prefix)
ski_name = simulation.prefix()

# The path to the simulation log file
remote_log_file_path = fs.join(simulation.remote_output_path, ski_name + "_log.txt")

# Check whether the log file exists
if not remote.is_file(remote_log_file_path): raise RuntimeError("The log file does not exist remotely")

# Read the log file
lines = remote.read_text_file(remote_log_file_path)

# Print the lines of the log file
for line in lines: print(line)

# -----------------------------------------------------------------
