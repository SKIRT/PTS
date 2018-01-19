#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

# Images

# -----------------------------------------------------------------

definition.add_flag("weighed", "plot weighed residuals", None)
definition.add_flag("distributions", "plot the residual distributions", False)
definition.add_flag("absolute", "plot absolute residuals", False)

# -----------------------------------------------------------------

# Extra flags
definition.add_flag("normalize", "normalize the images")
definition.add_flag("share_scale", "share the scales of the images")
definition.add_optional("scale_reference", "string", "name of the image to determine the scale for to use for the other images")
definition.add_flag("same_residuals_scale", "use the same scale for the residuals as for the observation and models")

# -----------------------------------------------------------------

definition.add_flag("write", "write out the processed frames, masks and regions", False)
definition.add_flag("show", "show the plot (default is automatic)", None)

# -----------------------------------------------------------------

# Add coordinates?
definition.add_flag("coordinates", "show the coordinates", True)

# -----------------------------------------------------------------
