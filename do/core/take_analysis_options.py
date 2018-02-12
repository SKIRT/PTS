#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.take_analysis_options Take the analysis options from a set of simulation files and set them to simulations with corresponding names.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.tools import filesystem as fs
from pts.core.basics.log import log
from pts.core.remote.host import find_host_ids
from pts.core.simulation.simulation import RemoteSimulation
from pts.core.simulation.remote import get_simulations_for_host
from pts.core.simulation.shower import compare_analysis

# -----------------------------------------------------------------

host_ids = find_host_ids()

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Remote host
definition.add_required("remote", "string", "remote host for which to find corresponding simulations", choices=host_ids)

# Show
definition.add_flag("show", "show the comparison between the analysis options")

# Read the command line arguments
config = parse_arguments("take_analysis_options", definition, description="Take the analysis options from a set of simulation files and set them to simulations with corresponding names")

# -----------------------------------------------------------------

#simulations = []
simulations = dict()

# Loop over the simulation files in the working directory
for path in fs.files_in_cwd(extension="sim"):

    # Load simulation
    simulation = RemoteSimulation.from_file(path)

    # Add the simulation to the list
    #simulations.append(simulation)

    # Add the simulation
    simulations[simulation.name] = simulation

# -----------------------------------------------------------------

# Get the simulation names
simulation_names = simulations.keys()

# -----------------------------------------------------------------

# Find corresponding simulations
simulations_host = get_simulations_for_host(config.remote, names=simulation_names, as_dict=True)

# -----------------------------------------------------------------

# Loop over the simulations
for name in simulations_host:

    # Debugging
    log.debug("Taking analysis options for simulation '" + name + "' ...")

    # Get simulation and reference simulation
    simulation = simulations_host[name]
    reference_simulation = simulations[name]

    # Show
    if config.show: compare_analysis(simulation, reference_simulation)

    # Set analysis options
    simulation.analysis = reference_simulation.analysis

    # Save the simulation
    simulation.save()

# -----------------------------------------------------------------
