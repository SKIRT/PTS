#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.remote.host import find_host_ids

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition(log_path="log", config_path="config")

# Add required arguments
definition.add_required("image", "string", "name of the image for which to run the preparation")

# Add optional arguments
definition.add_optional("exclude_filters", "string_list", "exclude the data for these filters from the procedure that brings all data to the same resolution and pixelscale")
definition.add_flag("steps", "write the results of intermediate steps")
definition.add_flag("visualise", "make visualisations")

# Remote preparation
definition.add_optional("remote", "string", "remote host on which to run the preparation", choices=find_host_ids())

# -----------------------------------------------------------------
