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

# Add optional
definition.add_optional("remote", str, "update on a remote system")
definition.add_optional("pubkey_password", str, "password of public key on the system")

# -----------------------------------------------------------------

