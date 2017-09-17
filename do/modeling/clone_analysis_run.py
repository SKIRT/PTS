#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.show_status Show the status page.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.tools import time
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

environment = load_modeling_environment_cwd()
analysis_runs = environment.analysis_runs

# -----------------------------------------------------------------

definition = ConfigurationDefinition()
definition.add_required("analysis_run", "string", "analysis run to clone", choices=analysis_runs.names)
config = parse_arguments("clone_analysis_run", definition)

# -----------------------------------------------------------------

# Get analysis run path
path = analysis_runs.get_path(config.analysis_run)

# -----------------------------------------------------------------

# Generate new analysis run name
new_name = time.unique_name()
new_path = fs.create_directory_in(environment.analysis_path, new_name)

# -----------------------------------------------------------------

clear_directories = ["out", "misc", "plot", "residuals", "heating", "extr", "colours", "attenuation"]
fs.copy_from_directory(path, new_path, exact_not_name=clear_directories, not_extension="sh")
fs.create_directories_in(new_path, clear_directories)

# -----------------------------------------------------------------
