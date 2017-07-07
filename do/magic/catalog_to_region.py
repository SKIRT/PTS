#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.catalog_to_regions Create a regions file from a catalog file.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import tables
from pts.magic.basics.coordinate import SkyCoordinate
from pts.magic.core.frame import Frame
from pts.magic.tools import statistics
from pts.core.tools import filesystem as fs
from pts.magic.catalog.extended import ExtendedSourceCatalog
from pts.magic.catalog.point import PointSourceCatalog
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()
definition.add_required("type", "string", "type of catalog (extended/point)", choices=["extended", "point"])
definition.add_required("path", "file_path", "name of the region file")
definition.add_required("frame", "file_path", "name of an image file")

#parser.add_argument("image", type=str, help="the name of the image file for which to create the region")
#parser.add_argument("fwhm", type=float, help="the FWHM of the stars (in pixels)")
#parser.add_argument("sigma_level", type=float, nargs='?', help="the sigma level", default=3.0)
#parser.add_argument("--color", type=str, help="the color", default="blue")

config = parse_arguments("catalog_to_regions", definition)

# -----------------------------------------------------------------

# Load the catalog
#catalog_name = os.path.splitext(os.path.basename(arguments.catalog))[0]
#catalog = tables.from_file(arguments.catalog)

# Open the frame
#frame = Frame.from_file(arguments.image)

# -----------------------------------------------------------------

catalog_name = fs.strip_extension(fs.name(config.path))

# -----------------------------------------------------------------

# Load the catalog
if config.type == "extended": catalog = ExtendedSourceCatalog.from_file(config.path)
elif config.type == "point": catalog = PointSourceCatalog.from_file(config.path)
else: raise ValueError("Invalid catalog type")

# -----------------------------------------------------------------



# -----------------------------------------------------------------