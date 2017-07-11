#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.modeling.config.generate_page import definition

# -----------------------------------------------------------------

# Set the modeling path
modeling_path = fs.cwd()

# -----------------------------------------------------------------

# Flags
definition.add_flag("use_session", "use remote python session to create the images", False)
definition.add_flag("show", "show the page", False)

# -----------------------------------------------------------------
