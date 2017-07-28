#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

# Create definition
definition = ConfigurationDefinition()

# Add optional
definition.add_positional_optional("script", "string", "name of the PTS do script for which to determine the dependencies")

# Add flags
definition.add_flag("modules", "show the PTS modules which import a given package", letter="m")
definition.add_flag("standard", "show import packages from the python standard library", letter="s")
definition.add_flag("version", "show the version numbers of the required packages", letter="v")
definition.add_flag("canopy", "show whether the package is available through Canopy")
definition.add_flag("pip", "show whether the package is available through pip")
definition.add_flag("conda", "show whether the package is available through the conda package manager")
definition.add_flag("description", "show a description")
definition.add_flag("subprojects", "show the dependencies per subproject", letter="j")

definition.add_flag("show", "show", True)
definition.add_flag("write", "write", True)

# -----------------------------------------------------------------
