#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.commands Contains the ModelingCommands class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

class ModelingCommands(list):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ModelingCommands, self).__init__(*args)

        # Set the path
        self.path = kwargs.pop("path", None)

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, filepath):

        """
        This function ...
        :param filepath:
        :return:
        """

        # Load the lines
        return cls(fs.get_lines(filepath), path=filepath)

    # -----------------------------------------------------------------

    def saveto(self, path, update_path=True):

        """
        This function ...
        :param path:
        :param update_path:
        :return:
        """

        # Write the lines
        fs.write_lines(path, self)

        # Update the path
        if update_path: self.path = path

    # -----------------------------------------------------------------

    def save(self):

        """
        This function ...
        :return:
        """

        if self.path is None: raise ValueError("The path is undefined")

        # Save
        self.saveto(self.path)

# -----------------------------------------------------------------
