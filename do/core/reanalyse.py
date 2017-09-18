#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.reanalyse Re-analyse a certain simulation.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.remote.host import all_host_ids
from pts.core.simulation.remote import get_simulation_for_host
from pts.core.launch.basicanalyser import steps, extraction, plotting, misc
from pts.core.launch.analyser import SimulationAnalyser

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()
definition.add_required("remote", "string", "remote host to mount", choices=all_host_ids())
definition.add_required("id", "positive_integer", "simulation ID")
definition.add_positional_optional("steps", "string_list", "re-analyse only certain steps", choices=steps, default=steps)

# Read the command line arguments
config = parse_arguments("reanalyse", definition, description="Re-analyse a certain simulation")

# -----------------------------------------------------------------

# Load the simulation
simulation = get_simulation_for_host(config.remote, config.id)

# Check whether retrieved
if not simulation.retrieved: raise ValueError("The simulation has not been retrieved yet")

# -----------------------------------------------------------------

# Flag indicating whether this simulation has been analysed or not
simulation.analysed = False

# Reset extraction
if extraction in config.steps: simulation.analysed_extraction = []

# Reset plotting
if plotting in config.steps: simulation.analysed_plotting = []

# Reset misc
if misc in config.steps: simulation.analysed_misc = []

# -----------------------------------------------------------------

# Create simulation analyser
analyser = SimulationAnalyser()

# Run the analyser on the simulation
analyser.run(simulation=simulation)

# -----------------------------------------------------------------
