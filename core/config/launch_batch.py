#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.host import find_host_ids
from pts.core.config.simulation.logging import definition as logging_definition
from pts.core.config.simulation.launch import definition
from pts.core.config.simulation.analysis import definition as analysis_definition

# -----------------------------------------------------------------

# Remote
definition.add_optional("remotes", "string_list", "remote host IDs to use", choices=find_host_ids(), default=find_host_ids())
definition.add_optional("extra_remote", "string", "remote host ID to use for the extra simulations", choices=find_host_ids())

definition.add_flag("shared_input", "whether the different simulations share their input folder", False)

definition.add_optional("cores_per_process", "integer", "the number of cores to use per process (if the parallelization is not specified by the user for a certain remote host)")

definition.add_flag("data_parallel", "enable data parallelization mode (if the parallelization is not specified by the user of the batch launcher for a certain remote host)")

definition.add_flag("group_simulations", "group multiple simulations in one job", False)

definition.add_optional("group_walltime", "real", "preferred walltime per job of grouped simulations")

definition.add_optional("timing_table_path", "file_path", "path to the timing table")
definition.add_optional("memory_table_path", "file_path", "path to the memory table")

# Logging
definition.import_section("logging", "logging options", logging_definition)

# Analysis options
definition.import_section("analysis", "simulation analysis options", analysis_definition)

# -----------------------------------------------------------------
