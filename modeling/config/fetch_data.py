#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import Configuration

# -----------------------------------------------------------------

# Create the configuration
config = Configuration(log_path="log")

# Add settings
config.add_section("database")
config.sections["database"].add_optional("username", str, "the username")
config.sections["database"].add_optional("password", str, "the password")

# -----------------------------------------------------------------
