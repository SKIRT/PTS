#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.list_data This ...

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import sys

# Import astronomical modules
from telarchive import archive_search

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition(log_path="log")

# The galaxy name
definition.add_required("galaxy", "string", "galaxy name")

# Get configuration
config = parse_arguments("list_data", definition)

# -----------------------------------------------------------------

archive_search.main(sys.argv)

# -----------------------------------------------------------------
