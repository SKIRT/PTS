#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

sigma_levels = [1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5.]

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")
definition.add_flag("show", "show the page", False)

# Sigma levels
definition.add_positional_optional("sigma_levels", "real_list", "different sigma levels for which to generate significance masks", sigma_levels)

# Flags
definition.add_flag("replot", "replot already existing figures", False) # default False because prep data doesn't really change (normally)

# ADVANCED
definition.add_optional("nopen_files", "positive_integer", "number of open files necessary to make the script work", 1024)

# -----------------------------------------------------------------
