#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.generation_reanalysis Reanalyse simulations of a generation.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.core.launch.basicanalyser import steps, extraction, plotting, misc
from pts.core.launch.analyser import SimulationAnalyser
from pts.core.launch.options import extraction_names, plotting_names, misc_names
from pts.core.tools import sequences

# -----------------------------------------------------------------

environment = load_modeling_environment_cwd()
runs = environment.fitting_runs

# -----------------------------------------------------------------

batch = "batch"
all_steps = steps + [batch]

# -----------------------------------------------------------------

# Create definition
definition = ConfigurationDefinition()

# The fitting run for which to explore the parameter space
if runs.empty: raise RuntimeError("No fitting runs are present")
elif runs.has_single: definition.add_fixed("name", "name of the fitting run", runs.single_name)
else: definition.add_required("name", "string", "name of the fitting run", choices=runs.names)

# Generations to remove
definition.add_required("generation", "string", "generation name")

# Reanalyse which steps?
definition.add_positional_optional("steps", "string_list", "re-analyse only certain steps", choices=all_steps, default=all_steps)
definition.add_positional_optional("features", "string_list", "re-analyse only certain features (if a single step is defined)")

# Create the configuration
config = parse_arguments("generation_reanalysis", definition)

# -----------------------------------------------------------------

# Load the fitting run and the generation
fitting_run = runs.load(config.name)
generation = fitting_run.get_generation(config.generation)

# -----------------------------------------------------------------

# Loop over the simulations
for simulation in generation.simulations:

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
        #if scaling in config.steps: simulation.analysed_scaling = False

    # Create simulation analyser
    analyser = SimulationAnalyser()

    # Run the analyser on the simulation
    analyser.run(simulation=simulation)

# -----------------------------------------------------------------
