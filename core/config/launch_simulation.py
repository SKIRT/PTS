#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.config.launch import definition
from pts.core.remote.host import find_host_ids

# -----------------------------------------------------------------

# Add required arguments
definition.add_required("ski", "file_path", "name/path of the ski file")

# Simulation settings
definition.add_flag("relative", "treats the given input and output paths as being relative to the ski/fski file")

# Flags
definition.add_flag("create_output", "create the output directory if it is not yet present")
definition.add_flag("show_progress", "show progress (for local simulation or remote simulation in attached mode)", True)

# -----------------------------------------------------------------

# Add positional arguments
definition.add_positional_optional("remote", "string", "remote host on which to run the simulation (if none is specified, the simulation is run locally", choices=find_host_ids())
definition.add_optional("cluster_name", "string", "name of the cluster", letter="c")

# Remote and parallelization
definition.add_optional("parallel", "integer_pair", "parallelization scheme (processes, threads)", letter="p")
definition.add_optional("walltime", "duration", "an estimate for the walltime of the simulation for the specified parallelization scheme")
definition.add_flag("data_parallel_local", "enable data parallelization for local execution", False)
definition.add_flag("data_parallel_remote", "set data parallelization for remote execution (None means automatic)", None)

# Check parallelization? (ADVANCED)
definition.add_flag("check_parallelization", "check the specified parallelization scheme", True)

# -----------------------------------------------------------------

definition.add_optional("finish_at", "string", "finish when this is encountered")
definition.add_optional("finish_after", "string", "finish after this is encountered")

# -----------------------------------------------------------------

definition.add_flag("debug_output", "output lines of SKIRT process", False)

# -----------------------------------------------------------------
