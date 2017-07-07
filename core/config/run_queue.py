#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add required
definition.add_required("name", "string", "name of the simulation queue")
definition.add_required("walltime", "duration", "walltime for the jobs")
definition.add_positional_optional("runid", "integer", "start with this run ID", default=0)

# -----------------------------------------------------------------
