#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

# Configuration
definition = ConfigurationDefinition()

# Galaxy name
definition.add_required("galaxy_name", "string", "the name of the galaxy")
definition.add_required("band", "string", "the band (GALEX or SDSS u/g/r/i/z)")

# Optional
definition.add_optional("remote", "string", "the remote host name", None)

# -----------------------------------------------------------------
