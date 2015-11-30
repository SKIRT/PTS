#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

"""
This module ...
"""

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import inspect
from distutils.spawn import find_executable

# -----------------------------------------------------------------

# The path to the root PTS directory
pts_root_dir = inspect.getfile(inspect.currentframe()).split("/pts")[0]

# The path to the PTS package directory (PTS/pts)
pts_package_dir = os.path.join(pts_root_dir, "pts")

# The path to the PTS user directory (PTS/user)
pts_user_dir = os.path.join(pts_root_dir, "user")

# -----------------------------------------------------------------

# The path to the root SKIRT directory
skirt_root_dir = find_executable("skirt").split("/release")[0]

# The path to the SKIRT repository
skirt_repo_dir = os.path.join(skirt_root_dir, "git")

# The path to the SKIRT release directory
skirt_release_dir = os.path.join(skirt_root_dir, "release")

# The path to the SKIRT run directory
skirt_run_dir = os.path.join(skirt_root_dir, "run")

# -----------------------------------------------------------------

def dependencies(module):

    """
    This function ...
    :param path:
    :return:
    """



# -----------------------------------------------------------------