#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.build.definition Contains the ModelDefinition class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

class ModelDefinition(object):
    
    """
    This class...
    """

    def __init__(self, name, path):

        """
        This function ..
        :param name:
        :param path:
        """

        # Set model name and directory path
        self.name = name
        self.path = path

        # Subdirectories
        self.stellar_path = fs.create_directory_in(self.path, "stellar")
        self.dust_path = fs.create_directory_in(self.path, "dust")

        # Other input
        self.input_path = fs.create_directory_in(self.path, "input")

    # -----------------------------------------------------------------

    @property
    def stellar_map_paths(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.stellar_path, recursive=True, exact_name="map", extension="fits")

    # -----------------------------------------------------------------

    @property
    def dust_map_paths(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.dust_path, recursive=True, exact_name="map", extension="fits")

    # -----------------------------------------------------------------

    @property
    def input_paths(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.input_path) + self.stellar_map_paths + self.dust_map_paths

# -----------------------------------------------------------------
