#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.remote.host import find_host_ids
from pts.core.config.launch import definition

# -----------------------------------------------------------------

# Flags
definition.add_flag("recursive", "look for ski files recursively")

# Settings for the remote
if len(find_host_ids()) > 0: definition.add_optional("remotes", "string_list", "remote host IDs to use", choices=find_host_ids(), default=find_host_ids())
else: definition.add_fixed("remotes", "remote hosts", [])

# Parallelization options
definition.add_optional("nnodes", "integer", "number of computing nodes to be used (for remote hosts that use a scheduling system, for other remotes the current load of the system will be probed)")

# Advanced options
definition.add_flag("shared_input", "whether the different simulations share their input folder", False)
definition.add_flag("group_simulations", "group multiple simulations in one job", False)
definition.add_flag("use_pts", "use PTS one the remote end to dynamicaly pick simulations for the job")
definition.add_optional("group_walltime", "real", "preferred walltime per job of grouped simulations")
definition.add_flag("progress_bar", "use progress bars to show progress")

# The timing and memory table
definition.add_optional("timing_table_path", "file_path", "path to the timing table")
definition.add_optional("memory_table_path", "file_path", "path to the memory table")

# Flags
definition.add_flag("add_timing_local", "add timing information for locally launched simulations", True)
definition.add_flag("add_memory_local", "add memory information for locally launched simulations", True)

# The analyser classes
definition.add_optional("analysers", "string_list", "analyser classes for the simulations")

definition.add_flag("relative_analysis_paths", "treat the analysis paths (extraction, plotting, misc) as bein relative to the respective ski files")

# -----------------------------------------------------------------
