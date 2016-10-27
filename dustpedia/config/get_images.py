#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.dustpedia.data.images import instruments

# -----------------------------------------------------------------

# Configuration
definition = ConfigurationDefinition()

# Galaxy name
definition.add_required("galaxy_name", "string", "the name of the galaxy")

# Database
definition.add_section("database", "DustPedia database credentials")
definition.add_optional("username", "string", "username")
definition.add_optional("password", "string", "password")

# Other options
definition.add_positional_optional("instruments", "string_list", "instruments for which to download the images", choices=instruments, default=instruments)

# -----------------------------------------------------------------
