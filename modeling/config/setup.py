#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.remote.host import find_host_ids

# -----------------------------------------------------------------

types = dict()
types["galaxy"] = "create a 3D model of a galaxy based on image data"
types["other"] = "create a model of an other object based on SED data"

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

# Add required settings
definition.add_required("type", "string", "type of modeling", choices=types)
definition.add_required("name", "string", "name given to the object")
definition.add_required("fitting_host_ids", "string_list", "remote hosts to use for performing simulations as part of the fitting", choices=find_host_ids())

# -----------------------------------------------------------------
