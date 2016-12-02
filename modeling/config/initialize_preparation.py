#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.basics.host import find_host_ids

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition(log_path="log", config_path="config")

# Add optional arguments
definition.add_optional("image", "string", "the name of the image for which to run the initialization")

# Add flags
definition.add_flag("visualise", "make visualisations")

# Remote source detection
definition.add_optional("remote", "string", "remote host on which to run the source finder", choices=find_host_ids())

# -----------------------------------------------------------------
