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
definition.add_positional_optional("host_id", "string", "install SKIRT/PTS remotely")
definition.add_optional("repository", "string", "repository name from which to clone (only possible when installing remotely and SKIRT/PTS is already installed locally)")

# Add flags
definition.add_flag("private", "use the private SKIRT/PTS repository")

# Advanced
definition.add_flag("force", "force re-installation when already present", letter="f")
definition.add_flag("force_conda", "force installation of conda and creation of a conda environment even when it is already present")

# Add flag
#definition.add_flag("all_remotes", "update on all remote hosts")

# For PTS:
definition.add_optional("python_name", "string", "name for the python environment for PTS", "python_pts")
definition.add_optional("python_version", "string", "version number for python", "2.7")

# -----------------------------------------------------------------
