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

# The path of the resulting mosaic image
definition.add_optional("band_id", "string", "ID of the band")
definition.add_optional("image_path", "file_path", "path to the resulting mosaic image file")

# Task id
definition.add_optional("host_id", "string", "host ID of the task", choices=find_host_ids())
definition.add_optional("task_id", "integer", "task ID")

# -----------------------------------------------------------------
