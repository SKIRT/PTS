#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.remote.host import find_host_ids

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()

# Add optional
definition.add_positional_optional("ski", "file_path", "ski file to be used (if not specified, all ski files in the current directory will be used)")

# Add flag
definition.add_flag("recursive", "look for ski files recursively")

# Add optional
definition.add_optional("remote", "string", "the remote host (no schedulers) on which to launch the simulations", choices=find_host_ids(schedulers=False))

# Add optional
definition.add_optional("nprocesses", "integer", "number of processes to use", 1)
definition.add_flag("data_parallel", "use data parallelization mode")

# -----------------------------------------------------------------
