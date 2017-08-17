#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.convolve Uniformize a set of images.

# -----------------------------------------------------------------

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.magic.core.list import NamedFrameList
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

# Create the definition
definition = ConfigurationDefinition()

# Select files
definition.add_optional("contains", "string_list", "only files whose name contain one of these strings")

# Unit
definition.add_optional("unit", "photometric_unit", "unit to which to convert all images")

# Flags
definition.add_flag("backup", "make backups", True)

# Parse the command line arguments
config = parse_arguments("uniformize", definition)

# -----------------------------------------------------------------

# Load frame list
frames = NamedFrameList.from_directory(config.path, contains=config.contains)

# -----------------------------------------------------------------

# Uniformize
frames.uniformize(unit=config.unit)

# -----------------------------------------------------------------

# Backup original files
if config.backup: fs.backup_files_in_directory(config.path, extension="fits", contains=config.contains)

# Save new frames
frames.write_to_directory(config.path, replace=True)

# -----------------------------------------------------------------
