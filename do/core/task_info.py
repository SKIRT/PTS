#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.task_info Show information about a certain remote PTS task.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.remote.host import find_host_ids

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add required
definition.add_required("remote", "string", "ID of the remote host", choices=find_host_ids())
definition.add_required("task_id", "integer", "ID of the PTS task")

# -----------------------------------------------------------------

# Parse the arguments into a configuration
config = parse_arguments("task_info", definition, description="Show information about a certain remote PTS task")

# -----------------------------------------------------------------

print("Not implemented yet")

# -----------------------------------------------------------------
