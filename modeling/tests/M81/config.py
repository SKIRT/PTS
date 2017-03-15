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
definition = ConfigurationDefinition(write_config=False)

# Optional settings
definition.add_optional("npackages", "positive_integer", "number of photon packages per wavelength for the reference simulation", int(1e6))

# For remote execution of reference simulation
definition.add_optional("host_ids", "string_list", "remote hosts to use for heavy computations (in order of preference)", choices=find_host_ids(schedulers=False))

# -----------------------------------------------------------------
