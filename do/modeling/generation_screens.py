#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.generation_screens View the screens of a certain generation.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.core.remote.mounter import mount_remote
from pts.core.tools import filesystem as fs
from pts.core.remote.remote import get_home_path, Remote
from pts.core.tools import formatting as fmt

# -----------------------------------------------------------------

environment = load_modeling_environment_cwd()
runs = environment.fitting_runs

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# The fitting run for which to explore the parameter space
if runs.empty: raise RuntimeError("No fitting runs are present")
elif runs.has_single: definition.add_fixed("name", "name of the fitting run", runs.single_name)
else: definition.add_required("name", "string", "name of the fitting run", choices=runs.names)

# Generations to remove
definition.add_required("generation", "string", "generation name")

# Flags
definition.add_flag("open_output", "open the screen output path")
definition.add_flag("show_output", "show the screen output")

# Get configuration
config = parse_arguments("generation_screens", definition, "View the screens of a certain generation")

# -----------------------------------------------------------------

# Load the fitting run
fitting_run = runs.load(config.name)

# Get the generation
generation = fitting_run.get_generation(config.generation)

# -----------------------------------------------------------------

# Check
if not generation.has_screen_scripts: raise RuntimeError("No screen sessions for this generation")
if not generation.has_assignment_table: raise RuntimeError("No assignment for this generation")

# -----------------------------------------------------------------

# Get number of simulations
nsimulations = generation.nsimulations

# -----------------------------------------------------------------

# Check
nscreens = generation.nscreens
if nscreens > 1 and config.open_output: raise ValueError("Multiple screens")

# -----------------------------------------------------------------

# Loop over the screen names
screen_filename = "screenlog.0"
for name in generation.screen_names:

    print("")
    print(fmt.green + fmt.underlined + name + fmt.reset)
    print("")
    script = generation.get_screen_script(name)
    print(script.output_path)

    # Show output
    #home_path = get_home_path(script.host_id)
    remote = Remote(host_id=script.host_id)
    home_path = remote.home_directory
    if config.open_output:

        # Mount
        mount_path = mount_remote(script.host_id)

        # Determine screen output path
        relative_output_path = fs.relative_to(script.output_path, home_path)
        output_path = fs.join(mount_path, relative_output_path)
        fs.open_directory(output_path)

    # Show screen output
    if config.show_output:

        script_filepath = fs.join(script.output_path, screen_filename)
        #for line in remote.read_lines(script_filepath):
        #    print(line)

        #lines = remote.read_first_lines(script_filepath, 100)
        #for line in lines: print(line)

        nsimulations = remote.get_noccurences(script_filepath, "Constructing a simulation from ski file")
        print(nsimulations)

# -----------------------------------------------------------------
