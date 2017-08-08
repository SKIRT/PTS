#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.developer.remove_pyc Remove all compiled python (.pyc) files.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import introspection
from pts.core.tools import filesystem as fs
from pts.core.basics.log import log

# -----------------------------------------------------------------

# Loop over all pyc files
for path in fs.files_in_path(introspection.pts_package_dir, recursive=True, extension="pyc"):

    # Inform the user
    log.info("Removing '" + path + "' ...")

    # Remove the file
    fs.remove_file(path)

# -----------------------------------------------------------------
