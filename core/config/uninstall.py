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

# Create definition
definition = ConfigurationDefinition()

# Add optional settings
definition.add_optional("remotes", "string_list", "remote host(s) on which to uninstall", default=find_host_ids(), choices=find_host_ids())
definition.add_positional_optional("skirt_and_or_pts", "string_tuple", "SKIRT and/or PTS", default=["skirt", "pts"], choices=["skirt", "pts"])

# Add flags
definition.add_flag("conda", "also remove conda installation")
definition.add_flag("qt", "also remove Qt installation")

# -----------------------------------------------------------------
