#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# The filter for decomposition
# TODO: let the choices be the filters for which there is data !!
definition.add_optional("filter", "filter", "filter for which to use the data for decomposition", default="IRAC I1")

# The method
definition.add_optional("method", "string", "method for decomposition", default="s4G", choices=["s4g", "fit", "imfit"])

# Add optional
choices = dict()
choices["azimuth"] = "azimuth angle and y flattening"
choices["tilt"] = "tilt angle and z flattening"
definition.add_optional("bulge_deprojection_method", "string", "method of deprojecting a 2D bulge with position angle difference w.r.t. the disk", choices=choices, default="tilt")

# -----------------------------------------------------------------
