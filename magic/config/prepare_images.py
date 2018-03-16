#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.magic.config.extract import definition as extraction_definition
from pts.magic.config.subtract_sky import definition as subtraction_definition
from pts.core.tools.parallelization import ncores

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

# The dataset or image
definition.add_positional_optional("dataset", "file_path", "name of the dataset file or image file")

# -----------------------------------------------------------------
