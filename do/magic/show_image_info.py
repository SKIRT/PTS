#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.show_image_info Show basic info of a FITS image.

# -----------------------------------------------------------------

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.magic.tools.info import get_image_info_strings_from_header_file
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

# Create the definition
definition = ConfigurationDefinition()
definition.add_required("file_path", "file_path", "name of the input image file")

# Parse the command line arguments
config = parse_arguments("show_image_info", definition)

# -----------------------------------------------------------------

filename = fs.strip_extension(fs.name(config.file_path))
for line in get_image_info_strings_from_header_file(filename, config.file_path, name=False): print(line)

# -----------------------------------------------------------------
