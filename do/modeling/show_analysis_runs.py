#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.show_analysis_runs Show the analysis runs.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.core.tools import formatting as fmt
from pts.core.tools.stringify import tostr
from pts.modeling.config.parameters import parameter_descriptions

# -----------------------------------------------------------------

# Determine the modeling path
environment = load_modeling_environment_cwd()
#runs = environment.analysis_runs
context = environment.analysis_context

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add flags
definition.add_flag("info", "show the analysis run info", False)
definition.add_optional("parameters", "string_list", "show the values of these parameters", choices=parameter_descriptions)

# Get configuration
config = parse_arguments("show_model", definition)

# -----------------------------------------------------------------

# Cannot define parameters and enable info
if config.info and config.parameters is not None: raise ValueError("Cannot specify parameters and enable info")

# -----------------------------------------------------------------

ignore_properties = ["name", "path"]

# -----------------------------------------------------------------

def show_basic_info(run):

    """
    This function ...
    :param run:
    :return:
    """

    print(run.info.to_string(line_prefix="    ", ignore_none=True, ignore=ignore_properties))

    # From config
    print("     - " + fmt.bold + "old scale heights: " + fmt.reset + tostr(run.config.old_scale_heights))

# -----------------------------------------------------------------

def show_wavelength_grid_info(run):

    """
    This function ...
    :param run:
    :return:
    """

    print("     - " + fmt.bold + "wavelength grid:" + fmt.reset)
    print("        - " + fmt.bold + "npoints: " + fmt.reset + tostr(run.nwavelengths))
    print("        - " + fmt.bold + "emission lines: " + fmt.reset + tostr(run.config.wg.add_emission_lines))

# -----------------------------------------------------------------

def show_dust_grid_info(run):

    """
    This function ...
    :param run:
    :return:
    """

    # From config
    print("     - " + fmt.bold + "dust grid:" + fmt.reset)
    print("        - " + fmt.bold + "type: " + fmt.reset + run.config.dg.grid_type)
    print("        - " + fmt.bold + "relative scale: " + fmt.reset + tostr(run.config.dg.scale))
    print("        - " + fmt.bold + "scale heights: " + fmt.reset + tostr(run.config.dg.scale_heights))
    if run.config.dg.grid_type == "bintree": print("        - " + fmt.bold + "min level: " + fmt.reset + tostr(run.config.dg.bintree_min_level))
    elif run.config.dg.grid_type == "octtree": print("        - " + fmt.bold + "min level: " + fmt.reset + tostr(run.config.dg.octtree_min_level))
    else: raise ValueError("Invalid grid type: " + run.config.dg.grid_type)
    print("        - " + fmt.bold + "maximum mass fraction: " + fmt.reset + tostr(run.config.dg.max_mass_fraction))

    # From dust grid object
    print("        - " + fmt.bold + "sample count: " + fmt.reset + tostr(run.dust_grid.sample_count))
    print("        - " + fmt.bold + "min x: " + fmt.reset + tostr(run.dust_grid.min_x))
    print("        - " + fmt.bold + "max x: " + fmt.reset + tostr(run.dust_grid.max_x))
    print("        - " + fmt.bold + "min y: " + fmt.reset + tostr(run.dust_grid.min_y))
    print("        - " + fmt.bold + "max y: " + fmt.reset + tostr(run.dust_grid.max_y))
    print("        - " + fmt.bold + "min z: " + fmt.reset + tostr(run.dust_grid.min_z))
    print("        - " + fmt.bold + "max z: " + fmt.reset + tostr(run.dust_grid.max_z))
    print("        - " + fmt.bold + "direction method: " + fmt.reset + tostr(run.dust_grid.direction_method))
    print("        - " + fmt.bold + "maximum optical depth: " + fmt.reset + tostr(run.dust_grid.max_optical_depth))
    print("        - " + fmt.bold + "maximum density dispersion fraction: " + fmt.reset + tostr(run.dust_grid.max_dens_disp_fraction))
    print("        - " + fmt.bold + "search method: " + fmt.reset + tostr(run.dust_grid.search_method))

# -----------------------------------------------------------------

def show_run_info(run):

    """
    This function ...
    :param run:
    :return:
    """

    # Show basic info
    show_basic_info(run)

    print("")

    # Show wavelength grid info
    show_wavelength_grid_info(run)

    print("")

    # Show dust grid info
    show_dust_grid_info(run)

# -----------------------------------------------------------------

# Loop over the cache host IDS
for host_id in context.cache_host_ids:

    print("")
    print(fmt.yellow + host_id.upper() + ":" + fmt.reset)
    print("")

    # Get the run names
    run_names = context.get_run_names_for_host_id(host_id)

    # Loop over the runs
    for name in run_names:

        # Show the name
        print(" - " + fmt.underlined + fmt.blue + name + fmt.reset)

        # Show the info
        if config.info:

            # Load the run
            run = context.get_cached_run(name)

            # Show the info
            print("")
            show_run_info(run)
            print("")

        # Show the parameters
        elif config.parameters is not None:

            # Load the run
            run = context.get_cached_run(name)

            print("")
            # Show the parameter values
            for name in config.parameters:
                print("   - " + fmt.bold + name + ": " + fmt.reset + tostr(run.info.parameter_values[name]))
            print("")

# -----------------------------------------------------------------

# Empty line to separate
if not (config.info or config.parameters is not None): print("")

# Show the local analysis runs
print(fmt.yellow + "LOCAL:" + fmt.reset)
print("")

for name in context.analysis_run_names:

    # Show the name
    print(" - " + fmt.underlined + fmt.blue + name + fmt.reset)

    # Show the info
    if config.info:

        # Load the run
        run = context.get_run(name)

        # Show the info
        print("")
        show_run_info(run)
        print("")

    # Show the parameters
    elif config.parameters is not None:

        # Load the run
        run = context.get_run(name)

        print("")
        # Show the parameter values
        for name in config.parameters:
            print("   - " + fmt.bold + name + ": " + fmt.reset + tostr(run.info.parameter_values[name]))
        print("")

# End with empty line
if not (config.info or config.parameters is not None): print("")

# -----------------------------------------------------------------
