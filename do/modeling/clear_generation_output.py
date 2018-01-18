#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.clear_generation_output Clear the output of all simulations of a certain generation.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import verify_modeling_cwd
from pts.modeling.fitting.run import FittingRuns
from pts.core.tools import filesystem as fs
from pts.core.basics.log import log

# -----------------------------------------------------------------

# Get the modeling commands
modeling_path = verify_modeling_cwd()
runs = FittingRuns(modeling_path)

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# The fitting run for which to explore the parameter space
if runs.empty: raise RuntimeError("No fitting runs are present")
elif runs.has_single: definition.add_fixed("name", "name of the fitting run", runs.single_name)
else: definition.add_required("name", "string", "name of the fitting run", choices=runs.names)

# Flags
definition.add_flag("backup", "make backups of non-empty directories", False)
definition.add_flag("adapt_simulations", "unset retrieved and analysed flags of the corresponding simulations", False)

# Generations to remove
definition.add_required("generation", "string", "generation to remove (none means all)")

# Get configuration
config = parse_arguments("clear_generation_output", definition)

# -----------------------------------------------------------------

# Load the fitting run
fitting_run = runs.load(config.name)

# Check
if not fitting_run.is_generation(config.generation): raise ValueError("Generation doesn't exist")

# Get generation path
generation_path = fitting_run.get_generation_path(config.generation)

# -----------------------------------------------------------------

# Loop over the simulation paths
for path, name in fs.directories_in_path(generation_path, returns=["path", "name"]):

    # Determine subpaths
    out_path = fs.join(path, "out")
    extr_path = fs.join(path, "extr")
    plot_path = fs.join(path, "plot")
    misc_path = fs.join(path, "misc")

    # Non-empty output
    if fs.is_directory(out_path) and not fs.is_empty(out_path):

        # Debugging
        log.debug("Clearing output directory of '" + name + "' simulation ...")
        if config.backup: fs.backup_directory(out_path)
        fs.clear_directory(out_path)

    # Non-empty extracted data
    if fs.is_directory(extr_path) and not fs.is_empty(extr_path):

        # Debugging
        log.debug("Clearing extraction directory of '" + name + "' simulation ...")
        if config.backup: fs.backup_directory(extr_path)
        fs.clear_directory(extr_path)

    # Non-empty plotting output
    if fs.is_directory(plot_path) and not fs.is_empty(plot_path):

        # Debugging
        log.debug("Clearing plotting directory of '" + name + "' simulation ...")
        if config.backup: fs.backup_directory(plot_path)
        fs.clear_directory(plot_path)

    # Non-empty misc output
    if fs.is_directory(misc_path) and not fs.is_empty(misc_path):

        # Debugging
        log.debug("Clearing miscellaneous output directory of '" + name + "' simulation ...")
        if config.backup: fs.backup_directory(misc_path)
        fs.clear_directory(misc_path)

# -----------------------------------------------------------------

# Clear chi-squared table
chi_squared_path = fitting_run.chi_squared_table_path_for_generation(config.generation)
chi_squared = fitting_run.chi_squared_table_for_generation(config.generation)

# Non-empty
if len(chi_squared) > 0:

    log.debug("Clearing chi-squared table ...")
    if config.backup: fs.backup_file(chi_squared_path)
    chi_squared.remove_all_rows()
    chi_squared.save()

# -----------------------------------------------------------------

# Adapt simulations?
if config.adapt_simulations:

    # Get the simulations for the generation
    generation = fitting_run.get_generation(config.generation)

    # Check if has assignment table
    if generation.has_assignment_table:

        # Get the simulations
        for simulation in generation.simulations:

            # Set flag
            changed = False

            # Unset retrieved
            changed = simulation.set_retrieved(False)

            # Unset analysed
            #print(simulation.retrieved, simulation.analysed, simulation.analysed_extraction, simulation.analysed_plotting, simulation.analysed_misc, simulation.analysed_batch, simulation.analysed_scaling, simulation.analysed_extra)
            changed |= simulation.set_analysed(False, all=True)

            # Has changed
            if changed:
                log.debug("Adapting retrieved and analysed flags for simulation '" + simulation.name + "' ...")
                simulation.save()

    # Give warning
    else: log.warning("Assignment table not found")

# -----------------------------------------------------------------
