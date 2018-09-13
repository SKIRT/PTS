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

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.core.tools import filesystem as fs
from pts.core.tools import introspection
from pts.core.tools import terminal
from pts.core.tools import formatting as fmt
from pts.core.tools import strings

# -----------------------------------------------------------------

# Load modeling environment
environment = load_modeling_environment_cwd()
runs = environment.analysis_runs

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

# The analysis run
if runs.empty: raise RuntimeError("No analysis runs present (yet)")
elif runs.has_single: definition.add_fixed("run", "name of the analysis run", runs.single_name)
else: definition.add_positional_optional("run", "string", "name of the analysis run", runs.last_name, runs.names)

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

# Show
for name, path in fs.directories_in_path(correlation_plots_path, returns=["name", "path"], sort=lambda name: int(name[0])):

    # Show
    print(fmt.bold + name + fmt.reset_bold)
    print("")

    # Determine path to the file containing the stilts command
    filepath = fs.get_filepath(path, "stilts.txt")

    # Get text
    #command = fs.get_text(filepath)
    #for line in lines: print(line)
    #print(command.replace("\\\n", ""))

    # Get the original command
    lines = fs.get_lines(filepath)
    command = ""
    for line in lines: command += line.split(" \\")[0].strip() + " "

    data_filepaths = strings.get_substrings(command, "/Users", ".dat", only_shortest=True)

    #print(command)
    #print(data_filepaths)
    for data_filepath in data_filepaths:

        rel_filepath = fs.relative_to(data_filepath, correlations_path)
        print(rel_filepath)

    #print(list(command))
    #exit()

    print("")

    # Execute the plotting command
    #terminal.execute(command)

# -----------------------------------------------------------------
