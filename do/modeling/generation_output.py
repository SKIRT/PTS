#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.generation_output Show the output of a generation.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.core.tools import filesystem as fs
from pts.core.tools import formatting as fmt

# -----------------------------------------------------------------

environment = load_modeling_environment_cwd()
runs = environment.fitting_runs

# -----------------------------------------------------------------

# Create definition
definition = ConfigurationDefinition()

# The fitting run for which to explore the parameter space
if runs.empty: raise RuntimeError("No fitting runs are present")
elif runs.has_single: definition.add_fixed("name", "name of the fitting run", runs.single_name)
else: definition.add_required("name", "string", "name of the fitting run", choices=runs.names)

# Generations to remove
definition.add_required("generation", "string", "generation name")

# Create the configuration
config = parse_arguments("generation_output", definition, "Show the output of a generation")

# -----------------------------------------------------------------

# Load the fitting run and the generation
fitting_run = runs.load(config.name)
generation = fitting_run.get_generation(config.generation)

# -----------------------------------------------------------------

print("")
print(fmt.yellow + "FILES" + fmt.reset)
print("")

# Show files
for filename in fs.files_in_path(generation.path, returns="name", extensions=True):
    print(" - " + fmt.bold + filename + fmt.reset)
print("")

# -----------------------------------------------------------------

print(fmt.yellow + "SIMULATIONS" + fmt.reset)
print("")

# -----------------------------------------------------------------

# Loop over the simulation names
for simulation_name in generation.simulation_names:

    # Check whether directory is present
    if not fs.contains_directory(generation.path, simulation_name):
        print(" - " + fmt.red + simulation_name + ": directory missing")
    else:

        simulation_path = fs.join(generation.path, simulation_name)
        ski_path = fs.join(simulation_path, environment.galaxy_name + ".ski")
        out_path = fs.join(simulation_path, "out")
        extr_path = fs.join(simulation_path, "extr")
        plot_path = fs.join(simulation_path, "plot")
        misc_path = fs.join(simulation_path, "misc")

        # Check presence
        has_ski = fs.is_file(ski_path)
        has_out = fs.is_directory(out_path)
        has_extr = fs.is_directory(extr_path)
        has_plot = fs.is_directory(plot_path)
        has_misc = fs.is_directory(misc_path)

        missing = []
        if not has_ski: missing.append("ski file")
        if not has_out: missing.append("output directory")
        if not has_extr: missing.append("extraction directory")
        if not has_plot: missing.append("plotting directory")
        if not has_misc: missing.append("misc directory")
        nmissing = len(missing)

        if nmissing == 0: print(" - " + fmt.green + simulation_name)
        else: print(" - " + fmt.red + simulation_name + ": " + ",".join(missing) + " missing")
print("")

# -----------------------------------------------------------------

# Define necessary output files
logfiles_name = "logfiles"
seds_name = "seds"
datacubes_name = "datacubes"

# -----------------------------------------------------------------

print(fmt.yellow + "OUTPUT" + fmt.reset)
print("")

# Initialize
missing_output = OrderedDict()
missing_output[logfiles_name] = []
missing_output[seds_name] = []
if generation.use_images: missing_output[datacubes_name] = []

# Loop over the simulations
for simulation_name in generation.simulation_names:

    # Get output
    output = generation.get_simulation_output(simulation_name)

    # Check missing files
    if not output.has_logfiles: missing_output[logfiles_name].append(simulation_name)
    if not output.has_seds: missing_output[seds_name].append(simulation_name)
    if generation.use_images and not output.has_total_images: missing_output[datacubes_name].append(simulation_name)

# -----------------------------------------------------------------

# Show
for output_type in missing_output:

    simulation_names = missing_output[output_type]
    nsimulations = len(simulation_names)

    if nsimulations == 0: print(" - " + fmt.green + output_type + fmt.reset)
    else: print(" - " + fmt.red + output_type + ": missing for " + ", ".join(simulation_names) + fmt.reset)
print("")

# -----------------------------------------------------------------

# Define necessary output files
timeline_name = "timeline"
memory_name = "memory"

# -----------------------------------------------------------------

print(fmt.yellow + "EXTRACTION" + fmt.reset)
print("")

# Initialize
missing_extraction = OrderedDict()
missing_extraction[timeline_name] = []
missing_extraction[memory_name] = []

# Loop over the simulations
for simulation_name in generation.simulation_names:

    # Get extraction output
    extraction = generation.get_extraction_output(simulation_name)

    # Check missing files
    if not extraction.has_timeline: missing_extraction[timeline_name].append(simulation_name)
    if not extraction.has_memory: missing_extraction[memory_name].append(simulation_name)

# -----------------------------------------------------------------

# Show
for output_type in missing_extraction:

    simulation_names = missing_extraction[output_type]
    nsimulations = len(simulation_names)

    if nsimulations == 0: print(" - " + fmt.green + output_type + fmt.reset)
    else: print(" - " + fmt.red + output_type + ": missing for " + ", ".join(simulation_names) + fmt.reset)
print("")

# -----------------------------------------------------------------

print(fmt.yellow + "PLOTTING" + fmt.reset)
print("")

# Initialize
missing_plotting = OrderedDict()
missing_plotting[seds_name] = []

# Loop over the simulations
for simulation_name in generation.simulation_names:

    # Get plotting output
    plotting = generation.get_plotting_output(simulation_name)

    # Check missing files
    if not plotting.has_seds: missing_plotting[seds_name].append(simulation_name)

# -----------------------------------------------------------------

# Show
for output_type in missing_plotting:

    simulation_names = missing_plotting[output_type]
    nsimulations = len(simulation_names)

    if nsimulations == 0: print(" - " + fmt.green + output_type + fmt.reset)
    else: print(" - " + fmt.red + output_type + ": missing for " + ", ".join(simulation_names) + fmt.reset)
print("")

# -----------------------------------------------------------------

# Define necessary output files and directories
differences_name = "differences"
fluxes_name = "fluxes"
images_name = "images"

# -----------------------------------------------------------------

print(fmt.yellow + "MISC" + fmt.reset)
print("")

# Initialize
missing_misc = OrderedDict()
missing_misc[differences_name] = []
missing_misc[fluxes_name] = []
if generation.use_images: missing_misc[images_name] = []

# Loop over the simulations
for simulation_name in generation.simulation_names:

    # Get misc output
    misc = generation.get_misc_output(simulation_name)

    # Check missing files and directories
    if generation.use_images:
        if not misc.has_image_fluxes: missing_misc[fluxes_name].append(simulation_name)
    else:
        if not misc.has_fluxes: missing_misc[fluxes_name].append(simulation_name)
    if generation.use_images and not misc.has_images_for_fluxes: missing_misc[images_name].append(simulation_name)
    if not misc.has_other_filename("differences.dat"): missing_misc[differences_name].append(simulation_name)

# -----------------------------------------------------------------

# Show
for output_type in missing_misc:

    simulation_names = missing_misc[output_type]
    nsimulations = len(simulation_names)

    if nsimulations == 0: print(" - " + fmt.green + output_type + fmt.reset)
    else: print(" - " + fmt.red + output_type + ": missing for " + ", ".join(simulation_names) + fmt.reset)
print("")

# -----------------------------------------------------------------
