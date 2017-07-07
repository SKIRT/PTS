#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.scale_regions Scale regions with a certain factor.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.tools import filesystem as fs
from pts.magic.region.list import load_region_list

# -----------------------------------------------------------------

definition = ConfigurationDefinition()
definition.add_required("path", "file_path", "region path list")
definition.add_required("factor", "positive_real", "scale factor")
config = parse_arguments("scale_regions", definition)

# -----------------------------------------------------------------

regions = load_region_list(config.path)

# -----------------------------------------------------------------

regions *= config.factor

# -----------------------------------------------------------------

# Remove the original file
fs.remove_file(config.path)

# Save the regions again
regions.saveto(config.path)

# -----------------------------------------------------------------
