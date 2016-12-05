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
definition.add_positional_optional("remote", "string", "install SKIRT/PTS remotely")
definition.add_optional("repository", "string", "repository name from which to clone (only possible when installing remotely and SKIRT/PTS is already installed locally)")

# Add flags
definition.add_flag("private", "use the private SKIRT/PTS repository")
definition.add_flag("force", "force re-installation when already present", letter="f")

# -----------------------------------------------------------------
