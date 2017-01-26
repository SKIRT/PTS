#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.remote.host import find_host_ids
from pts.core.config.simulation.logging import definition as logging_definition
from pts.core.config.simulation.launch import definition
from pts.core.config.simulation.analysis import definition as analysis_definition

# -----------------------------------------------------------------

# Flags
definition.add_flag("recursive", "look for ski files recursively")

# Settings for the remote
definition.add_optional("remotes", "string_list", "remote host IDs to use", choices=find_host_ids(), default=find_host_ids())

# Parallelization options
definition.add_optional("nnodes", "integer", "number of computing nodes to be used (for remote hosts that use a scheduling system, for other remotes the current load of the system will be probed)")

# Advanced options
definition.add_flag("shared_input", "whether the different simulations share their input folder", False)
definition.add_flag("group_simulations", "group multiple simulations in one job", False)
definition.add_flag("use_pts", "use PTS one the remote end to dynamicaly pick simulations for the job")
definition.add_optional("group_walltime", "real", "preferred walltime per job of grouped simulations")

# The timing and memory table
definition.add_optional("timing_table_path", "file_path", "path to the timing table")
definition.add_optional("memory_table_path", "file_path", "path to the memory table")

# Logging options
definition.import_section("logging", "logging options", logging_definition)

# Analysis options
definition.import_section("analysis", "simulation analysis options", analysis_definition)

# The analyser classes
definition.add_optional("analysers", "string_list", "analyser classes for the simulations")

# -----------------------------------------------------------------
