#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.remote.host import find_host_ids

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

# Remote
definition.add_optional("remote", "string", "remote host on which to execute the calculations", choices=find_host_ids())

# Method
definition.add_optional("method", "string", "fitting method", "grid", choices=["grid", "genetic"])

# -----------------------------------------------------------------
