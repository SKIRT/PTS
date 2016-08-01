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

# Add settings
definition.add_section("database", "options for the DustPedia database connection")
definition.sections["database"].add_optional("username", "string", "the username")
definition.sections["database"].add_optional("password", "string", "the password")

# -----------------------------------------------------------------
