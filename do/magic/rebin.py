#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.rebin Rebin a set of images to the same pixel grid.

# -----------------------------------------------------------------

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.magic.core.list import NamedFrameList
from pts.core.tools import filesystem as fs
from pts.magic.basics.coordinatesystem import CoordinateSystem

# -----------------------------------------------------------------

# Create the definition
definition = ConfigurationDefinition()

# Select files
definition.add_optional("contains", "string_list", "only files whose name contain one of these strings")

# Options
definition.add_optional("name", "string", "rebin to this name")
definition.add_optional("wcs", "file_path", "path to the file from which to take the WCS")

# Flags
definition.add_flag("backup", "make backups", True)

# Parse the command line arguments
config = parse_arguments("rebin", definition)

# -----------------------------------------------------------------

# Load frame list
frames = NamedFrameList.from_directory(config.path, contains=config.contains)

# -----------------------------------------------------------------

# Rebin
if config.name is not None: frames.rebin_to_name(config.name)
elif config.wcs is not None:
    wcs = CoordinateSystem.from_file(config.wcs)
    frames.rebin_to_wcs(wcs)
else: frames.rebin_to_highest_pixelscale()

# -----------------------------------------------------------------

# Backup original files
if config.backup: fs.backup_files_in_directory(config.path, extension="fits", contains=config.contains)

# Save new frames
frames.write_to_directory(config.path, replace=True)

# -----------------------------------------------------------------
