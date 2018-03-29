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
from pts.modeling.analysis.heating.component import contributions

# -----------------------------------------------------------------

# Determine the modeling path
environment = load_modeling_environment_cwd()
context = environment.analysis_context

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add positional optional
definition.add_positional_optional("runs", "string_list", "names of analysis runs")

# Add flags
definition.add_flag("info", "show the analysis run info", False)
definition.add_flag("cached", "include cached analysis runs", True)
definition.add_flag("all_info", "show all info", False)

# Different sections
definition.add_flag("basic", "show basic info", True)
definition.add_flag("model", "show model info", False)
definition.add_flag("wavelength", "show wavelength grid info", True)
definition.add_flag("grid", "show dust grid info", True)
definition.add_flag("launch", "show launch info", False)
definition.add_flag("parallelization", "show parallelization", False)
definition.add_flag("timing", "show timing info", False)
definition.add_flag("memory", "show memory usage info", False)
definition.add_flag("directory", "show directory info", False)
definition.add_flag("heating", "show heating launch info", False)
definition.add_flag("residuals", "show residuals analysis info", False)
definition.add_flag("colours", "show colours analysis info", False)
definition.add_flag("attenuation", "show attenuation analysis info", False)
definition.add_flag("maps", "show maps analysis info", False)

# Add optional
definition.add_optional("parameters", "string_list", "show the values of these parameters", choices=parameter_descriptions)

# Get configuration
config = parse_arguments("show_model", definition)

# -----------------------------------------------------------------

# Cannot define parameters and enable info
if config.info and config.parameters is not None: raise ValueError("Cannot specify parameters and enable info")

# -----------------------------------------------------------------

ignore_properties = ["name", "path"]

# -----------------------------------------------------------------

# All info
if config.all_info:
    if not config.info: raise ValueError("All info is enabled but showing info is disabled")
    config.basic = True
    config.model = True
    config.wavelength = True
    config.grid = True
    config.launch = True
    config.parallelization = True
    config.timing = True
    config.memory = True
    config.directory = True
    config.heating = True
    config.residuals = True
    config.colours = True
    config.attenuation = True
    config.maps = True

# -----------------------------------------------------------------

def show_basic_info(run):

    """
    This function ...
    :param run:
    :return:
    """

    # Show the contents of the info file
    print(run.info.to_string(line_prefix="    ", ignore_none=True, ignore=ignore_properties))

    # From config
    print("     - " + fmt.bold + "old scale heights: " + fmt.reset + tostr(run.config.old_scale_heights))

# -----------------------------------------------------------------

def show_model_info(run):

    """
    This function ...
    :param run:
    :return:
    """

    print("     - " + fmt.bold + "model:" + fmt.reset)
    print("")

    print("       - old stars:")
    print("         - map name: " + run.model_old_map_name)
    print("         - methods: " + tostr(run.old_map_methods))
    print("         - origins: " + tostr(run.old_map_origins))

    print("       - young stars:")
    print("         - map name: " + run.model_young_map_name)
    print("         - methods: " + tostr(run.young_map_methods))
    print("         - origins: " + tostr(run.young_map_origins))

    print("       - ionizing stars:")
    print("         - map name: " + run.model_ionizing_map_name)
    print("         - methods: " + tostr(run.ionizing_map_methods))
    print("         - origins: " + tostr(run.ionizing_map_origins))

    print("       - dust:")
    print("         - map name: " + run.model_dust_map_name)
    print("         - methods: " + tostr(run.dust_map_methods))
    print("         - origins: " + tostr(run.dust_map_origins))

# -----------------------------------------------------------------

def show_wavelength_grid_info(run):

    """
    This function ...
    :param run:
    :return:
    """

    print("     - " + fmt.bold + "wavelength grid:" + fmt.reset)
    print("        - " + fmt.bold + "number of points: " + fmt.reset + tostr(run.nwavelengths))
    print("        - " + fmt.bold + "emission lines: " + fmt.reset + tostr(run.config.wg.add_emission_lines))

# -----------------------------------------------------------------

def show_dust_grid_info(run):

    """
    This function ...
    :param run:
    :return:
    """

    # Title
    print("     - " + fmt.bold + "dust grid:" + fmt.reset)

    # From log file or tree file
    print("        - " + fmt.bold + "number of cells: " + fmt.reset + tostr(run.ncells))

    # From config
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

def show_launch_info(run):

    """
    This function ...
    :param run:
    :return:
    """

    # Get the ski file
    ski = run.ski_file

    # Show the launch options
    print("     - " + fmt.bold + "launch info:" + fmt.reset)
    print("        - " + fmt.bold + "number of photon packages: " + fmt.reset + tostr(ski.packages()))
    print("        - " + fmt.bold + "dust self-absorption: " + fmt.reset + tostr(ski.dustselfabsorption()))
    print("        - " + fmt.bold + "transient heating: " + fmt.reset + tostr(ski.transient_dust_emissivity))
    print("        - " + fmt.bold + "has output: " + fmt.reset + tostr(run.has_output))
    print("        - " + fmt.bold + "has extracted data: " + fmt.reset + tostr(run.has_extracted))
    print("        - " + fmt.bold + "has plots: " + fmt.reset + tostr(run.has_plots))
    print("        - " + fmt.bold + "has misc output: " + fmt.reset + tostr(run.has_misc))

# -----------------------------------------------------------------

def show_parallelization_info(run):

    """
    This function ...
    :param run:
    :return:
    """

    # Get the log file
    logfile = run.logfile

    # Show the parallelization info
    print("     - " + fmt.bold + "parallelization options:" + fmt.reset)
    print("        - " + fmt.bold + "number of processes: " + fmt.reset + tostr(logfile.nprocesses))
    print("        - " + fmt.bold + "number of threads: " + fmt.reset + tostr(logfile.nthreads))
    print("        - " + fmt.bold + "data parallelization: " + fmt.reset + tostr(logfile.data_parallel))

# -----------------------------------------------------------------

def show_timing_info(run):

    """
    Thisn function ...
    :param run:
    :return:
    """

    # Show the timing info
    print("     - " + fmt.bold + "timing:" + fmt.reset)
    print("        - " + fmt.bold + "total: " + fmt.reset + tostr(run.timeline.total))
    print("        - " + fmt.bold + "setup: " + fmt.reset + tostr(run.timeline.setup))
    print("        - " + fmt.bold + "stellar: " + fmt.reset + tostr(run.timeline.stellar))
    print("        - " + fmt.bold + "spectra: " + fmt.reset + tostr(run.timeline.spectra))
    print("        - " + fmt.bold + "dust: " + fmt.reset + tostr(run.timeline.dust))
    print("        - " + fmt.bold + "writing: " + fmt.reset + tostr(run.timeline.writing))
    print("        - " + fmt.bold + "communication: " + fmt.reset + tostr(run.timeline.communication))
    print("        - " + fmt.bold + "waiting: " + fmt.reset + tostr(run.timeline.waiting))

# -----------------------------------------------------------------

def show_memory_info(run):

    """
    This function ...
    :param run:
    :return:
    """

    # Show the memory info
    print("     - " + fmt.bold + "memory:" + fmt.reset)
    print("        - " + fmt.bold + "peak: " + fmt.reset + tostr(run.memory.peak))
    print("        - " + fmt.bold + "peak per process: " + fmt.reset + tostr(run.memory.peak_per_process))

# -----------------------------------------------------------------

def show_directory_info(run):

    """
    This function ...
    :param run:
    :return:
    """

    # Show the launch options
    print("     - " + fmt.bold + "directory info:" + fmt.reset)
    print("        - " + fmt.bold + "number of files: " + fmt.reset + tostr(run.nfiles))
    print("        - " + fmt.bold + "directory size: " + fmt.reset + tostr(run.disk_space))

# -----------------------------------------------------------------

def show_heating_info(run):

    """
    Thisf unction ...
    :param run:
    :return:
    """

    print("     - " + fmt.bold + "heating simulations:" + fmt.reset)
    print("")

    # BASIC
    print("        - " + fmt.bold + "basic:" + fmt.reset)
    print("           - " + fmt.bold + "old scale heights: " + fmt.reset + tostr(run.heating_config.old_scale_heights))

    print("")

    # Loop over the contributions
    npackages = None
    selfabsorption = None
    transient_heating = None
    for contribution in contributions:

        # Get the ski path
        #ski_path = run.heating_ski_path_for_contribution(contribution)
        #ski = SkiFile(ski_path)
        ski = run.get_ski_for_contribution(contribution)

        if npackages is None: npackages = ski.packages()
        elif ski.packages() != npackages: raise RuntimeError("")

        if selfabsorption is None: selfabsorption = ski.dustselfabsorption()
        elif ski.dustselfabsorption() != selfabsorption: raise RuntimeError("")

        if transient_heating is None: transient_heating = ski.transient_dust_emissivity
        elif ski.transient_dust_emissivity != transient_heating: raise RuntimeError("")

        # Get the output path
        #output_path = run.heating_output_path_for_contribution(contribution)

    # LAUNCH INFO
    print("        - " + fmt.bold + "launch info:" + fmt.reset)
    print("           - " + fmt.bold + "number of photon packages: " + fmt.reset + tostr(npackages))
    print("           - " + fmt.bold + "dust self-absorption: " + fmt.reset + tostr(selfabsorption))
    print("           - " + fmt.bold + "transient heating: " + fmt.reset + tostr(transient_heating))
    print("")

    # Loop over the contributions
    print("        - " + fmt.bold + "finished:" + fmt.reset)

    for contribution in contributions:

        output_path = run.output_path_for_contribution(contribution)
        print("           - " + fmt.bold + contribution + ": " + fmt.reset + tostr(output_path))

# -----------------------------------------------------------------

def show_residuals_info(run):

    """
    This function ...
    :param run:
    :return:
    """

    print("     - " + fmt.bold + "residuals:" + fmt.reset)
    print("")

    print("      - images: " + tostr(run.residual_image_names, delimiter=", "))

# -----------------------------------------------------------------

def show_colours_info(run):

    """
    This function ...
    :param run:
    :return:
    """

    print("     - " + fmt.bold + "colours:" + fmt.reset)
    print("")

    print("      - colours: " + tostr(run.colour_names, delimiter=", "))

# -----------------------------------------------------------------

def show_attenuation_info(run):

    """
    This function ...
    :param run:
    :return:
    """

    print("     - " + fmt.bold + "attenuation:" + fmt.reset)
    print("")

# -----------------------------------------------------------------

def show_maps_info(run):

    """
    Thisf unction ...
    :param run:
    :return:
    """

    print("     - " + fmt.bold + "maps:" + fmt.reset)
    print("")

    if run.has_maps_colours: print("      - colours: " + str(run.ncolour_maps) + " maps")
    if run.has_maps_ssfr: print("      - ssfr: " + str(run.nssfr_maps) + " maps")
    if run.has_maps_tir: print("      - tir: " + str(run.ntir_maps) + " maps")
    if run.has_maps_attenuation: print("      - attenuation: " + str(run.nattenuation_maps) + " maps")
    if run.has_maps_old: print("      - old: " + str(run.nold_maps) + " maps")
    if run.has_maps_dust: print("      - dust:" + str(run.dust_maps) + " maps")
    if run.has_maps_young: print("      - young: " + str(run.nyoung_maps) + " maps")
    if run.has_maps_ionizing: print("      - ionizing:" + str(run.nionizing_maps) + " maps")

# -----------------------------------------------------------------

def show_run_info(run):

    """
    This function ...
    :param run:
    :return:
    """

    # Show basic info
    if config.basic:
        show_basic_info(run)
        print("")

    # Show model info
    if config.model:
        show_model_info(run)
        print("")

    # Show wavelength grid info
    if config.wavelength:
        show_wavelength_grid_info(run)
        print("")

    # Show dust grid info
    if config.grid:
        show_dust_grid_info(run)
        print("")

    # Show launch info
    if config.launch:
        show_launch_info(run)
        print("")

    # Show parallelization
    if config.parallelization and run.has_logfile:
        show_parallelization_info(run)
        print("")

    # Show timing info
    if config.timing and run.has_timeline:
        show_timing_info(run)
        print("")

    # Show memory info
    if config.memory and run.has_memory:
        show_memory_info(run)
        print("")

    # Show directory info
    if config.directory:
        show_directory_info(run)
        print("")

    # Show heating launch info
    if config.heating and run.has_heating:
        show_heating_info(run)
        print("")

    # Show residuals analysis info
    if config.residuals and run.has_residuals:
        show_residuals_info(run)
        print("")

    # Show colours analysis info
    if config.colours and run.has_colours:
        show_colours_info(run)
        print("")

    # Show attenuation analysis info
    if config.attenuation and run.has_attenuation:
        show_attenuation_info(run)
        print("")

    # Show maps analysis info
    if config.maps and run.has_maps:
        show_maps_info(run)
        print("")

# -----------------------------------------------------------------

if config.cached:

    # Loop over the cache host IDS
    for host_id in context.cache_host_ids:

        print("")
        print(fmt.yellow + host_id.upper() + ":" + fmt.reset)
        print("")

        # Get the run names
        run_names = context.get_run_names_for_host_id(host_id)

        # Loop over the runs
        for name in run_names:

            # Check in runs
            if config.runs is not None and name not in config.runs: continue

            # Show the name
            print(" - " + fmt.underlined + fmt.blue + name + fmt.reset)

            # Show the info
            if config.info:

                # Load the run
                run = context.get_cached_run(name)

                # Show the info
                print("")
                show_run_info(run)
                #print("")

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

# Loop over the names
for name in context.analysis_run_names:

    # Check in runs
    if config.runs is not None and name not in config.runs: continue

    # Show the name
    print(" - " + fmt.underlined + fmt.blue + name + fmt.reset)

    # Show the info
    if config.info:

        # Load the run
        run = context.get_run(name)

        # Show the info
        print("")
        show_run_info(run)
        #print("")

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
