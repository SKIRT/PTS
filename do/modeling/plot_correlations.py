#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.plot_correlations Plot correlations of a galaxy analysis model.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.core.tools import filesystem as fs
from pts.core.tools import introspection
from pts.core.tools import terminal
from pts.core.tools import formatting as fmt
from pts.core.tools import strings
from pts.core.tools import sequences
from pts.core.basics.log import log
from pts.core.basics.plot import plotting_formats

# -----------------------------------------------------------------

# Load modeling environment
environment = load_modeling_environment_cwd()
runs = environment.analysis_runs

# -----------------------------------------------------------------

default_format = "png"

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

# The analysis run
if runs.empty: raise RuntimeError("No analysis runs present (yet)")
elif runs.has_single: definition.add_fixed("run", "name of the analysis run", runs.single_name)
else: definition.add_positional_optional("run", "string", "name of the analysis run", runs.last_name, runs.names)

# Show?
definition.add_flag("show", "show plot instead of writing")

# Plotting format
definition.add_optional("format", "string", "plotting format", default=default_format, choices=plotting_formats)

# Get configuration
config = parse_arguments("plot_correlations", definition)

# -----------------------------------------------------------------

# Load the analysis run
context = environment.analysis_context
analysis_run = context.get_run(config.run)
correlations_path = analysis_run.correlations_path

# -----------------------------------------------------------------

# Determine the path to the correlation plots directory
correlation_plots_path = fs.join(introspection.pts_dat_dir("modeling"), "CorrelationPlots")

# -----------------------------------------------------------------

print("")

# Loop over the different subdirectories of the correlations path
for correlation_name, correlation_path in fs.directories_in_path(correlations_path, returns=["name", "path"]):

    # Show correlation
    print(fmt.bold + "CORRELATION: " + correlation_name.upper() + fmt.reset_bold)
    print("")

    # Correlation plots path
    plots_reference_path = fs.join(correlation_plots_path, correlation_name)
    if not fs.is_directory(plots_reference_path): continue # no plotting routines yet for this correlation

    # Make plots path
    plots_path = fs.create_directory_in(correlation_path, "plots")

    # Show
    for reference_name, path in fs.directories_in_path(plots_reference_path, returns=["name", "path"], sort=lambda name: int(name[0])):

        # Show
        print("    " + fmt.bold + reference_name + fmt.reset_bold)
        print("")

        # Determine path to the file containing the stilts command
        filepath = fs.get_filepath(path, "stilts.txt")

        # Get the original command
        lines = fs.get_lines(filepath)
        command = ""
        for line in lines: command += line.split(" \\")[0].strip() + " "

        # Get filepaths
        data_filepaths = strings.get_substrings(command, "/Users", ".dat", only_shortest=True)

        # Split in internal and external
        correlation_data_filepaths = OrderedDict()
        external_data_filepaths = OrderedDict()
        for data_filepath in data_filepaths:
            if fs.contains_path(correlation_path, data_filepath):
                rel_filepath = fs.relative_to(data_filepath, correlation_path)
                correlation_data_filepaths[rel_filepath] = data_filepath
            else:
                rel_filepath = fs.relative_to(data_filepath, correlations_path)
                external_data_filepaths[rel_filepath] = data_filepath

        #print("    CORR", correlation_data_filepaths)
        #print("    EXT", external_data_filepaths)

        # Only external files (e.g. bd_ratio.dat in cartesian space): skip
        ncorrelation_files = len(correlation_data_filepaths)
        nexternal_files = len(external_data_filepaths)
        nfiles = ncorrelation_files + nexternal_files
        if nfiles == 0: raise IOError("Something is wrong with command in '" + filepath + "': cannot find data filepaths")
        if ncorrelation_files == 0:
            log.warning("No correlation data files in '" + filepath + "' command: skipping ...")
            continue

        # Create directory within plots path for this reference plot
        plot_path = fs.create_directory_in(plots_path, reference_name)

        mode = None
        the_data_filepath = None
        general_data_filepaths = None

        # At least one data file is a cells_xxx.dat one
        if strings.any_startswith(correlation_data_filepaths.keys(), "cells"):

            # Split in specific files and general files
            general_data_filepaths = strings.get_all_not_startswith(correlation_data_filepaths.keys(), "cells")
            cells_data_filepaths = strings.get_all_startswith(correlation_data_filepaths.keys(), "cells")

            # Check
            if sequences.all_equal(cells_data_filepaths): the_data_filepath = cells_data_filepaths[0]
            else: raise RuntimeError("Don't know what to do") # pairs of files are used for this plot (e.g. cells_a.dat and cells_b.dat)?

            # Set mode
            mode = "cells"

        # At least one data file is a pixels_xxx.dat one
        elif strings.any_startswith(correlation_data_filepaths.keys(), "pixels"):

            # Split in specific files and general files
            general_data_filepaths = strings.get_all_not_startswith(correlation_data_filepaths.keys(), "pixels")
            pixels_data_filepaths = strings.get_all_startswith(correlation_data_filepaths.keys(), "pixels")

            # Check
            if sequences.all_equal(pixels_data_filepaths): the_data_filepath = pixels_data_filepaths[0]
            else: raise RuntimeError("Don't know what to do") # pairs

            # Set mode
            mode = "pixels"

        # No cells or pixels, but clearly one data file per plot (except for external)
        elif sequences.all_equal(correlation_data_filepaths):

            the_data_filepath = correlation_data_filepaths.keys()[0]
            general_data_filepaths = []
            mode = "all"

        # We cannot know what to do
        else: raise RuntimeError("Don't know what to do")

        # Set full reference data filepath
        full_original_data_filepath = correlation_data_filepaths[the_data_filepath]

        # Set startswith
        startswith = mode if mode != "all" else None

        # Loop over the files
        for data_filename, data_filepath in fs.files_in_path(correlation_path, extension="dat", startswith=startswith, returns=["name", "path"]):

            #print(data_filename)
            if mode == "all": name = data_filename
            else:
                if data_filename.startswith(mode + "_"): name = data_filename.split(mode + "_")[1]
                else: name = data_filename

            # Debugging
            log.debug("Plotting " + name + " ...")

            # Determine plot filename
            plot_filename = name + "." + config.format
            if fs.has_file(plot_path, plot_filename): continue # Plot already exists

            # Create new command
            new_command = command.replace(full_original_data_filepath, data_filepath)
            if not config.show: new_command += " out='" + plot_filename + "'"

            #print(new_command)

            # Execute the plotting command
            terminal.execute(new_command, show_output=True, cwd=plot_path)

        #print(list(command))
        #exit()

        print("")

# -----------------------------------------------------------------
