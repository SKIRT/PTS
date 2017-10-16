#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.data Contains the SimulationData class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .output import SimulationOutput

# -----------------------------------------------------------------

class SimulationData(object):

    """
    This class ...
    """

    def __init__(self, output):

        """
        The constructor ...
        :param output:
        """

        # The simulation output
        self.output = output

    # -----------------------------------------------------------------

    @classmethod
    def from_directory(cls, path, prefix=None):

        """
        Thinsfunction ...
        :param path:
        :param prefix:
        :return:
        """

        # Get the output
        output = SimulationOutput.from_directory(path, prefix)

        # Create
        return cls(output)

    # -----------------------------------------------------------------

    @classmethod
    def from_cwd(cls, prefix=None):

        """
        This fucntion ...
        :param prefix:
        :return:
        """

        # Get the output
        output = SimulationOutput.from_cwd(prefix)

        # Create
        return cls(output)

    # -----------------------------------------------------------------

    @classmethod
    def from_paths(cls, paths):

        """
        This function ...
        :param paths:
        :return:
        """

        # Get the output
        output = SimulationOutput.from_paths(paths)

        # Create
        return cls(output)

    # -----------------------------------------------------------------

    def __str__(self):

        """
        This function ...
        :return:
        """

        return self.to_string()

    # -----------------------------------------------------------------

    def to_string(self, line_prefix=""):

        """
        This function ...
        :param line_prefix:
        :return:
        """

    # -----------------------------------------------------------------

    def show(self, line_prefix=""):

        """
        This function ...
        :param line_prefix:
        :return:
        """

        print(self.to_string(line_prefix=line_prefix))

# -----------------------------------------------------------------
