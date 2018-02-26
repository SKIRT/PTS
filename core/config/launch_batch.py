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

all_host_ids = find_host_ids()

# -----------------------------------------------------------------

# Flags
definition.add_flag("recursive", "look for ski files recursively")

# Load from queues
definition.add_flag("load_queues", "load queue files found in the working directory")

# Settings for the remote
if len(find_host_ids()) > 0: definition.add_optional("remotes", "string_list", "remote host IDs to use", choices=all_host_ids, default=all_host_ids)
else: definition.add_fixed("remotes", "remote hosts", [])

# Parallelization options
definition.add_optional("nnodes", "integer", "number of computing nodes to be used (for remote hosts that use a scheduling system, for other remotes the current load of the system will be probed)")
definition.add_optional("nnodes_per_host", "string_integer_dictionary", "number of computing nodes to be used, per scheduling remote")

# Advanced options
definition.add_flag("shared_input", "whether the different simulations share their input folder", None)
definition.add_flag("group_simulations", "group multiple simulations in one job", False)
definition.add_optional("group_walltime", "real", "preferred walltime per job of grouped simulations")
definition.add_flag("progress_bar", "use progress bars to show progress")

# The timing and memory table
definition.add_optional("timing_table_path", "file_path", "path to the timing table")
definition.add_optional("memory_table_path", "file_path", "path to the memory table")
definition.add_optional("runtimes_plot_path", "directory_path", "path for plotting runtimes (from runtime estimation)")
definition.add_flag("same_requirements", "set this flag when all simulations have approximately the same requirements", False)

# Flags
definition.add_flag("add_timing", True)
definition.add_flag("add_memory", True)
definition.add_flag("add_timing_local", "add timing information for locally launched simulations", True)
definition.add_flag("add_memory_local", "add memory information for locally launched simulations", True)

# The analyser classes
definition.add_optional("analysers", "string_list", "analyser classes for the simulations")

# -----------------------------------------------------------------

definition.add_flag("relative_analysis_paths", "treat the analysis paths (extraction, plotting, misc) as being relative to the respective ski files")

# -----------------------------------------------------------------

# Data parallel flags
definition.add_flag("data_parallel_local", "enable data parallelization for local execution", False)
definition.add_flag("data_parallel_remote", "set data parallelization for remote execution (None means automatic)", None)

# Check parallelization? (ADVANCED)
definition.add_flag("check_parallelization", "check the specified parallelization scheme", True)

# -----------------------------------------------------------------

# Write?
definition.add_flag("write", "enable writing", False)
definition.add_flag("write_assignment", "write assignment scheme", True)
definition.add_flag("write_queues", "write queues", True)

# -----------------------------------------------------------------

# Test mode
definition.add_flag("test", "test mode", False)

# -----------------------------------------------------------------

# VERY ADVANCED
definition.add_flag("all_sockets", "use all sockets, not just the determined number of 'free' sockets", False)
definition.add_optional("nsockets", "positive_integer", "use this number of sockets")
definition.add_flag("allow_multisocket_processes", "allow using multiple sockets per process", False)

# -----------------------------------------------------------------

definition.add_flag("cancel_scheduling_after_fail", "cancel the scheduling of simulations still in the queue of a remote host after the scheduling of a particular simulation was unsuccesful", True)
definition.add_flag("cancel_launching_after_fail", "cancel the launching of simulations still in the local queue after launching a particular simulation was unsuccesful", True)

# -----------------------------------------------------------------

# Showing
definition.add_flag("show_info", "show queue info", True)
definition.add_flag("show_status", "show remote status", True)
definition.add_flag("show_parallelizations", "show parallelization schemes", True)

# -----------------------------------------------------------------

# Other
definition.add_flag("clear_existing", "clear existing directories on the remote (otherwise error is thrown)", False)

# -----------------------------------------------------------------
