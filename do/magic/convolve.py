#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.convolve Convolve a set of images to the same resolution.

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

# Options
definition.add_optional("name", "string", "convolve to the resolution of the image with this name")
definition.add_optional("filter", "filter", "convolve to this filter's resolution")
definition.add_optional("fwhm", "angle", "specify the FWHM if the filter has a variable FWHM")

# Flags
definition.add_flag("backup", "make backups", True)

# Parse the command line arguments
config = parse_arguments("convolve", definition)

# -----------------------------------------------------------------

# Load frame list
frames = NamedFrameList.from_directory(config.path, contains=config.contains)

# -----------------------------------------------------------------

# Convolve
if config.name is not None: frames.convolve_to_name(config.name)
elif config.filter is not None: frames.convolve_to_filter(config.filter, fwhm=config.fwhm)
else: frames.convolve_to_highest_fwhm()

# -----------------------------------------------------------------

# Backup original files
if config.backup: fs.backup_files_in_directory(config.path, extension="fits", contains=config.contains)

# Save new frames
frames.write_to_directory(config.path, replace=True)

# -----------------------------------------------------------------
