#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.config.generate_page import definition
from pts.modeling.core.environment import verify_modeling_cwd

# -----------------------------------------------------------------

# Set the modeling path
modeling_path = verify_modeling_cwd()

# -----------------------------------------------------------------

# Flags
definition.add_flag("use_session", "use remote python session to create the images", False)

# Group
definition.add_flag("group_observatories", "group the images with observatories", False)

# Flags
definition.add_flag("thumbnails", "add map thumbnails", True)
definition.add_optional("thumbnail_height", "positive_integer", "height of the thumbnails (in pixels)", 50)
definition.add_flag("previews", "add previews of the maps when hovering over the thumbnails", True)

# -----------------------------------------------------------------
