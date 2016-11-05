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

# Required settings
#definition.add_required("galaxy_name", "string", "galaxy name")
definition.add_required("wcs", "file_path", "FITS file with the desired WCS")
definition.add_required("parameters", "directory_path", "path parameters directory")

# The output directory
definition.add_optional("output", "directory_path", "output directory", letter="o")

# -----------------------------------------------------------------
