#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.list_ukidss List UKIDSS images for a particular galaxy.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.magic.services.ukidss import UKIDSS

# -----------------------------------------------------------------

definition = ConfigurationDefinition()
definition.add_required("galaxy_name", "string", "galaxy name")
config = parse_arguments("list_ukidss", definition)

# -----------------------------------------------------------------

ukidss = UKIDSS()

#urls = ukidss.get_image_urls(config.galaxy_name)
#print(urls)

#names = ukidss.get_image_names(config.galaxy_name)
#print(names)

#filters = ukidss.get_image_filters(config.galaxy_name)
#print([str(fltr) for fltr in filters])

table = ukidss.get_image_table(config.galaxy_name)
print(table)

# -----------------------------------------------------------------
