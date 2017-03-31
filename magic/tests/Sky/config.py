#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.remote.host import find_host_ids
from pts.modeling.tests.base import possible_free_parameters, default_free_parameters
from pts.modeling.fitting.configuration import genetic_definition

# -----------------------------------------------------------------

# Create definition
definition = ConfigurationDefinition(write_config=False)

# -----------------------------------------------------------------
