#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

# Flags
definition.add_flag("interactive", "set the image paths manually (pass a list of image paths or ask for image paths interactively")
definition.add_flag("recursive", "search for image files recursively (within subdirectories)")

# Optional settings
definition.add_optional("error_suffix", "string", "the suffix used for the names of the error maps")
definition.add_optional("contains", "string", "only load files which contain this string")
definition.add_optional("not_contains", "string", "files with names that contain this string will be ignored")
definition.add_optional("exclude", "string", "exclude files with this (exact) name")

# -----------------------------------------------------------------
