#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

from pts.core.config.simulation.analysis import definition as analysis_definition
from pts.core.config.simulation.logging import definition as logging_definition
from pts.core.config.simulation.launch import definition
from pts.core.remote.host import find_host_ids

# -----------------------------------------------------------------

# Add required arguments
definition.add_required("ski", "file_path", "the name/path of the ski file")

# Simulation settings
definition.add_flag("relative", "treats the given input and output paths as being relative to the ski/fski file")

# -----------------------------------------------------------------

# Add positional arguments
definition.add_positional_optional("remote", "string", "the remote host on which to run the simulation (if none is specified, the simulation is run locally", choices=find_host_ids())

# Remote and parallelization
definition.add_optional("cluster", "string", "the name of the cluster", letter="c")
definition.add_optional("parallel", "integer_tuple", "the parallelization scheme (processes, threads)", letter="p")
definition.add_optional("walltime", "duration", "an estimate for the walltime of the simulation for the specified parallelization scheme")

# Logging options
definition.import_section("logging", "logging options", logging_definition)

# Analysis
definition.import_section("analysis", "analysis options", analysis_definition)

# -----------------------------------------------------------------
