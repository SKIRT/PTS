#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.environment Contains the ModelingEnvironment class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from .history import ModelingHistory

# -----------------------------------------------------------------

class ModelingEnvironment(object):

    """
    This class...
    """

    def __init__(self, modeling_path):

        """
        This function ...
        :param modeling_path
        """

        # Set the modeling base path
        self.path = modeling_path

        # Determine the path to the modeling configuration file
        self.config_file_path = fs.join(self.path, "modeling.cfg")

        # Determine the path to the modeling history file
        self.history_file_path = fs.join(self.path, "history.dat")

        # Initialize the history file
        if not fs.is_file(self.history_file_path):
            history = ModelingHistory()
            history.saveto(self.history_file_path)

        # Get the full paths to the necessary subdirectories and CREATE THEM
        self.fit_path = fs.create_directory_in(self.path, "fit")
        self.analysis_path = fs.create_directory_in(self.path, "analysis")
        self.reports_path = fs.create_directory_in(self.path, "reports")
        self.visualisation_path = fs.create_directory_in(self.path, "visualisation")
        self.plot_path = fs.create_directory_in(self.path, "plot")
        self.log_path = fs.create_directory_in(self.path, "log")
        self.config_path = fs.create_directory_in(self.path, "config")
        self.show_path = fs.create_directory_in(self.path, "show")
        self.build_path = fs.create_directory_in(self.path, "build")

# -----------------------------------------------------------------

class GalaxyModelingEnvironment(ModelingEnvironment):

    """
    This function ...
    """

    def __init__(self, modeling_path):

        """
        This function ...
        :param modeling_path:
        """

        # Call the constructor of the base class
        super(GalaxyModelingEnvironment, self).__init__(modeling_path)

# -----------------------------------------------------------------

class SEDModelingEnvironment(ModelingEnvironment):

    """
    This class ...
    """

    def __init__(self, modeling_path):

        """
        This function ...
        :param modeling_path:
        """

        # Call the constructor of the base class
        super(SEDModelingEnvironment, self).__init__(modeling_path)

# -----------------------------------------------------------------
