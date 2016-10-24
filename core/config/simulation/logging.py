#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

# Logging options
definition.add_flag("brief", "brief console logging", False, letter="b")
definition.add_flag("verbose", "verbose logging mode", False, letter="v")
definition.add_flag("memory", "memory logging", False, letter="m")
definition.add_flag("allocation", "memory (de)allocation logging", False, letter="a")
definition.add_optional("allocation_limit", "real", "memory (de)allocation logging lower limit (GB)", 1e-5, letter="l")
definition.add_flag("emulate", "emulate the simulation while limiting computation", letter="e")

# -----------------------------------------------------------------
