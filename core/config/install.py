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

definition.add_flag("finish", "finish previously initiated installation")

# Add flag
#definition.add_flag("all_remotes", "update on all remote hosts")

# For PTS:
definition.add_optional("python_name", "string_no_spaces", "name for the python environment for PTS (also the alias for the python executable)", "python_pts")
definition.add_optional("pip_name", "string_no_spaces", "name for the pip alias", "pip_pts")
definition.add_optional("jupyter_name", "string_no_spaces", "name for the jupyter executable", "jupyter_pts")
definition.add_optional("python_version", "string", "version number for python", "2.7")

# For SKIRT: re-install Qt
definition.add_optional("qt_version", "string", "version of Qt to install")
definition.add_flag("reinstall_qt", "re-install Qt")

# -----------------------------------------------------------------
