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
from pts.evolve.solve.extremizer import genetic_definition
from pts.core.tools.random import skirt_seed

# -----------------------------------------------------------------

# Create definition
definition = ConfigurationDefinition(write_config=False)

# Galaxy name
definition.add_optional("galaxy", "string", "galaxy name", "M77")

# Reference test image data
definition.add_optional("reference_path", "directory_path", "use the image data in this directory")
definition.add_optional("reference_test", "string", "use the image data of this previous test")

# -----------------------------------------------------------------
