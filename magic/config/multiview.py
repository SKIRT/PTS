#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.magic.config.view import definition

# -----------------------------------------------------------------

default_extensions = ["fits"]

# -----------------------------------------------------------------

# Image paths
definition.add_positional_optional("image_paths", "filepath_list", "image paths")
definition.add_optional("extensions", "string_list", "extensions of files to load", default_extensions)
definition.add_flag("recursive", "recursively load from directory", False)
definition.add_optional("contains", "string_list", "strings that have to be contained in the name of the files to be loaded")
definition.add_optional("not_contains", "string_list", "strings that cannot be contained in the name of the files to be loaded")
definition.add_optional("exact_name", "string_list", "exact name(s)")
definition.add_optional("exact_not_name", "string_list", "exact not name(s)")

# Regions
definition.add_optional("regions_prefix", "string", "prefix of regions filenames", "")
definition.add_optional("regions_suffix", "string", "suffix of regions filenames", "")
definition.add_optional("regions", "string", "exact filename for the regions files (or absolute path)")
definition.add_optional("regions_extension", "string", "extension of regions files", "reg")

# -----------------------------------------------------------------

definition.add_optional("width", "positive_integer", "image width", 300)
definition.add_optional("height", "positive_integer", "image height", 300)

# -----------------------------------------------------------------

# Preload
definition.add_flag("preload_all", "preload all images", False)
definition.add_optional("preload", "string_list", "names for which to preload the image")
definition.add_flag("dynamic", "create the viewers dynamically", False)

# -----------------------------------------------------------------

definition.add_flag("info", "add info about the images", False)

# -----------------------------------------------------------------

definition.add_optional("max_ncharacters_title", "positive_integer", "maximum number of characters in the titles before breaking line", 45)

# -----------------------------------------------------------------
