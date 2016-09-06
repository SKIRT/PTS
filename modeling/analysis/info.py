#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.info Contains the AnalysisRunInfo class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.configuration import load_mapping, write_mapping
from ...core.basics.map import Map

# -----------------------------------------------------------------

class AnalysisRunInfo(object):
    
    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        """

        self.name = None
        self.path = None
        self.generation_name = None
        self.simulation_name = None
        self.parameter_values = None
        self.chi_squared = None

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :return:
        """

        map = Map()
        with open(path, 'r') as infofile: load_mapping(infofile, map)

        # Set the info
        info = cls()
        info.name = map.name
        info.path = map.path
        info.generation_name = map.generation_name
        info.simulation_name = map.simulation_name
        info.parameter_values = map.parameter_values
        info.chi_squared = map.chi_squared

        # Return the info
        return info

    # -----------------------------------------------------------------

    def save(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        map = Map()
        map.name = self.name
        map.path = self.path
        map.generation_name = self.generation_name
        map.simulation_name = self.simulation_name
        map.parameter_valus = self.parameter_values
        map.chi_squared = self.chi_squared
        with open(path, 'w') as infofile: write_mapping(infofile, map)

# -----------------------------------------------------------------
