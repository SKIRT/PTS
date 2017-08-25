#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.remote.host import smb_host_ids

# -----------------------------------------------------------------

default_host_name = "www"

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# Upload
definition.add_optional("host_name", "string", "remote host name", default_host_name, choices=smb_host_ids())

# Flags
definition.add_flag("regenerate", "regenerate the pages", False)
definition.add_flag("replot", "make plots again", False)

# -----------------------------------------------------------------
