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

# Create the configuration definition
definition = ConfigurationDefinition()

# Add required
definition.add_required("remote", "string", "name of the remote host", choices=find_host_ids())
definition.add_positional_optional("id", "positive_integer", "simulation ID")
definition.add_optional("name", "string", "simulation name")

# Additional relative error
definition.add_optional("additional_error", "percentage", "additional percentual error for the observed flux points")

# Flags
definition.add_flag("write", "write the results", True)

# -----------------------------------------------------------------
