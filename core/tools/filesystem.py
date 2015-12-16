#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

"""
This module provides useful functions for manipulating the local file system.
"""

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os

# -----------------------------------------------------------------

def create_directory(path, recursive=False):

    """
    This function ...
    :param path:
    :param recursive:
    :return:
    """

    # Check whether the directory does not exist yet
    if not os.path.isdir(path):

        # Create the directory, recursively or not
        if recursive: os.makedirs(path)
        else: os.mkdir(path)

# -----------------------------------------------------------------

def create_directories(paths, recursive=False):
    
    """
    This function ...
    """
    
    # Loop over the different paths in the list
    for path in paths: create_directory(path)

# -----------------------------------------------------------------
