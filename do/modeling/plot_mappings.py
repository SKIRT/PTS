#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.plot_mappings Plot MAPPINGS examples.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.tools import filesystem as fs
from pts.modeling.misc.playground import plot_mappings_examples
from pts.core.basics.log import log

# -----------------------------------------------------------------

# Create configuration
definition = ConfigurationDefinition()
definition.add_positional_optional("directory", "directory_path", "output directory", fs.cwd())
config = parse_arguments("plot_mappings", definition)

# -----------------------------------------------------------------

# Create temp path
temp_path = fs.create_directory_in(config.directory, "seds")

# Plot
plot_mappings_examples(config.directory, temp_path=temp_path)

# -----------------------------------------------------------------
