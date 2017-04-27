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

# Galaxy
#definition.add_positional_optional("galaxy", "string", "galaxy to use images from")

# Number of frames
#definition.add_positional_optional("filters", "filter_list", "filters for which to use the images"

# -----------------------------------------------------------------
