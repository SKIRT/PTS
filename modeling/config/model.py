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

# Add flags
definition.add_flag("check_hosts", "check the availability of the remote hosts", True)
definition.add_flag("deploy", "deploy SKIRT and PTS where necessary", True)
definition.add_flag("check_versions", "check versions of SKIRT and PTS where necessary", True)

# Advanced settings
definition.add_flag("local", "keep computationaly heavy computations local")
definition.add_optional("remote", "string", "remote host for computationally heavy computations (overrule the modeling configuration)", choices=find_host_ids())
definition.add_flag("fitting_local", "launch the simulations as part of the fitting locally (overrule the modeling configuration")
definition.add_optional("fitting_remote", "string", "remote host for the fitting (overrule the modeling configuration)", choices=find_host_ids())

# -----------------------------------------------------------------
