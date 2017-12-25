#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.remote.host import smb_host_ids
from pts.modeling.config.component import definition

# -----------------------------------------------------------------

default_host_name = "www"

# -----------------------------------------------------------------

definition = definition.copy()

# Upload
definition.add_optional("host_name", "string", "remote host name", default_host_name, choices=smb_host_ids())

# Flags
definition.add_flag("generate", "generate the pages", True)
definition.add_flag("regenerate", "regenerate all pages", False)
definition.add_flag("replot", "make plots again", False)
definition.add_flag("show", "show after publishing", False)

definition.add_flag("regenerate_index", "regenerate the index page", False)

# Add detail pages
definition.add_flag("details", "add detailed pages", True)

# -----------------------------------------------------------------
