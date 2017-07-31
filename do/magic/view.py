#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.view View an image with JS9.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.magic.view.html import
from pts.core.tools import html
from pts.core.tools import filesystem as fs
from pts.core.tools import browser

# -----------------------------------------------------------------

definition = ConfigurationDefinition()
definition.add_required("image", "file_path", "image path")
definition.add_optional("regions", "file_path", "regions file path")
definition.add_optional("mask", "file_path", "mask file path")

# -----------------------------------------------------------------

browser.open_html()

# -----------------------------------------------------------------
