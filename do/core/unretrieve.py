#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.clear Clear retrieved, crashed, cancelled and aborted simulations

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import argparse

# Import the relevant PTS classes and modules
from pts.core.simulation.simulation import RemoteSimulation
from pts.core.tools import filesystem
from pts.core.basics.remote import Remote
from pts.core.tools import logging, time, parsing

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument("remote", nargs='?', default=None, help="the name of the remote host to connect to")
parser.add_argument("ids", type=parsing.simulation_ids, help="unretrieve the simulations with these ID's")
parser.add_argument("--keep", action="store_true", help="add this option to make sure the output is kept remotely after a subsequent retrieve attempt")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Create the remote execution environment
remote = Remote()
remote.setup(arguments.remote)



for path in filesystem.files_in_path(filesystem.cwd(), extension="sim"):
    
    # Create simulation
    sim = RemoteSimulation.from_file(path)
    
    if not sim.retrieved: continue
    
    local_output_path = sim.output_path
    
    remote_simulation_path = sim.remote_simulation_path
    remote_output_path = sim.remote_output_path
    remote_input_path = sim.remote_input_path
    
    if not remote.is_directory(remote_simulation_path): remote.create_directory(remote_simulation_path)
    if not remote.is_directory(remote_output_path): remote.create_directory(remote_output_path)
    
    #nancy.upload(local_output_path, remote_output_path)
    
    for filepath in filesystem.files_in_path(local_output_path):
        nancy.upload(filepath, remote_output_path, show_output=True)
    
    sim.retrieved = False
    
    if arguments.keep: sim.remove_remote_output = False
    
    sim.save()

# -----------------------------------------------------------------
