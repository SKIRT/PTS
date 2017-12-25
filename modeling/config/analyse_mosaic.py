#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.remote.host import find_host_ids
from pts.modeling.config.component import definition

# -----------------------------------------------------------------

definition = definition.copy()

# The path of the resulting mosaic image
definition.add_optional("band_id", "string", "ID of the band")
definition.add_optional("image_path", "file_path", "path to the resulting mosaic image file")
definition.add_optional("images_path", "directory_path", "path to a directory with mosaic images files")
definition.add_optional("out_path", "directory_path", "path of the out directory")

# Task id
definition.add_optional("host_id", "string", "host ID of the task", choices=find_host_ids())
definition.add_optional("task_id", "integer", "task ID")

# -----------------------------------------------------------------
