#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()

# Add optional
definition.add_positional_optional("host_id", "string", "update on a remote system")

# Add flags
definition.add_flag("dependencies", "also update the dependencies", False)

# Add flag
#definition.add_flag("all_remotes", "update on all remote hosts")

definition.add_flag("conda", "update conda", False)

# -----------------------------------------------------------------
