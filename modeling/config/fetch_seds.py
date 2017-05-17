#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.config.fetch import definition
from pts.magic.services.seds import catalog_names

# -----------------------------------------------------------------

definition.add_optional("catalogs", "string_list", "names of the catalogs to query", default=catalog_names, choices=catalog_names)

# -----------------------------------------------------------------
