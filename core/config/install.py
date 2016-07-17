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

# Add optional arguments
definition.add_optional("remote", str, "install SKIRT/PTS remotely")

# Add flags
definition.add_flag("private", "use the private SKIRT/PTS repository")

# -----------------------------------------------------------------
