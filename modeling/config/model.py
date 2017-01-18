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
definition = ConfigurationDefinition()

# Add flags
definition.add_flag("check_hosts", "check the availability of the remote hosts", True)
definition.add_flag("deploy", "deploy SKIRT and PTS where necessary", True)
definition.add_flag("check_versions", "check versions of SKIRT and PTS where necessary", True)

# -----------------------------------------------------------------
