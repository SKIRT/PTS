#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# Selections
definition.add_optional("old")
definition.add_optional("young")
definition.add_optional("ionizing")
definition.add_optional("dust")

# Anti-selections
definition.add_optional("not_old")
definition.add_optional("not_young")
definition.add_optional("not_ionizing")
definition.add_optional("not_dust")

# Flags
definition.add_flag("all_old")
definition.add_flag("all_young")
definition.add_flag("all_ionizing")
definition.add_flag("all_dust")
definition.add_flag("all")

# -----------------------------------------------------------------
