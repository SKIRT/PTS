#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.remote.host import find_host_ids
from pts.modeling.modeling.base import fitting_methods, default_fitting_method

# -----------------------------------------------------------------

types = dict()
types["galaxy"] = "Create a 3D model of a galaxy based on image data"
types["sed"] = "Model any object based on observed SED data"
types["images"] = "Model any object based on observational images and a ski file"

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(write_config=False)

# Add required settings
definition.add_required("type", "string", "type of modeling", choices=types)
definition.add_required("name", "string", "name given to the object")
definition.add_optional("fitting_method", "string", "fitting method", default_fitting_method, choices=fitting_methods)
definition.add_optional("fitting_host_ids", "string_list", "remote hosts to use for performing simulations as part of the fitting", choices=find_host_ids())

# -----------------------------------------------------------------
