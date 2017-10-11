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
from pts.core.launch.options import extraction_names, plotting_names, misc_names
from pts.core.tools import sequences

# -----------------------------------------------------------------

batch = "batch"
scaling = "scaling"

# -----------------------------------------------------------------

all_steps = steps + [batch, scaling]

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()
definition.add_required("remote", "string", "remote host to mount", choices=all_host_ids())
definition.add_required("id", "positive_integer", "simulation ID")
definition.add_positional_optional("steps", "string_list", "re-analyse only certain steps", choices=all_steps, default=all_steps)
definition.add_positional_optional("features", "string_list", "re-analyse only certain features (if a single step is defined)")

# Read the command line arguments
config = parse_arguments("reanalyse", definition, description="Re-analyse a certain simulation")

# -----------------------------------------------------------------

# Load the simulation
simulation = get_simulation_for_host(config.remote, config.id)

# Check whether retrieved
if not simulation.retrieved: raise ValueError("The simulation has not been retrieved yet")

# -----------------------------------------------------------------

# Only one step defined
if len(config.steps) == 1:

    # Flag indicating whether this simulation has been analysed or not
    simulation.analysed = False

    # Get the step
    step = config.steps[0]

    # Extraction
    if step == extraction:

        # All features
        if config.features is None: simulation.analysed_extraction = []

        # Select features
        else: simulation.analysed_extraction = sequences.elements_not_in_other(extraction_names, config.features, check_existing=True)

    # Plotting
    elif step == plotting:

        # All features
        if config.features is None: simulation.analysed_plotting = []

        # Select features
        else: simulation.analysed_plotting = sequences.elements_not_in_other(plotting_names, config.features, check_existing=True)

    # Misc
    elif step == misc:

        # All features
        if config.features is None: simulation.analysed_misc = []

        # Select features
        else: simulation.analysed_misc = sequences.elements_not_in_other(misc_names, config.features, check_existing=True)

    # Batch
    elif step == batch:

        # Check whether features are not defined
        if config.features is not None: raise ValueError("Cannot define features for the batch simulation analysis")

        # Set analysed_batch flag to False
        simulation.analysed_batch = False

    # Scaling
    elif step == scaling:

        # Check whether features are not defined
        if config.features is not None: raise ValueError("Cannot define features for scaling simulation analysis")

        # Set analysed_scaling flag to False
        simulation.analysed_scaling = False

# -----------------------------------------------------------------

# Multiple steps defined
else:

    # Check whether features are not defined
    if config.features is not None: raise ValueError("Features cannot be specified with multiple steps")

    # Flag indicating whether this simulation has been analysed or not
    simulation.analysed = False

    # Reset extraction
    if extraction in config.steps: simulation.analysed_extraction = []

    # Reset plotting
    if plotting in config.steps: simulation.analysed_plotting = []

    # Reset misc
    if misc in config.steps: simulation.analysed_misc = []

    # Reset batch
    if batch in config.steps: simulation.analysed_batch = False

    # Reset scaling
    if scaling in config.steps: simulation.analysed_scaling = False

# -----------------------------------------------------------------

# Create simulation analyser
analyser = SimulationAnalyser()

# Run the analyser on the simulation
analyser.run(simulation=simulation)

# -----------------------------------------------------------------
