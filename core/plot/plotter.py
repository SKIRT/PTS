#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.plot.plotter Contains the abstract Plotter class, from which all other plotter classes inherit.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta, abstractmethod, abstractproperty

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable
from ..basics.log import log

# -----------------------------------------------------------------

class Plotter(Configurable):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(Plotter, self).__init__(None)

        # -- Attributes --

        # The table containing the input data
        self.table = None

        # The data structure to contain the input information in plottable format
        self.data = []

    # -----------------------------------------------------------------

    @staticmethod
    @abstractproperty
    def default_input():

        """
        This function ...
        :return:
        """

        return ""

    # -----------------------------------------------------------------

    def run(self, table, output_path):

        """
        This function should be called to invoke the plotting routine ...
        :param table:
        :param output_path:
        :return:
        """

        # 1. Call the setup function
        self.setup(table, output_path)

        # 2. Prepare the input data into plottable format
        self.prepare_data()

        # 3. Create the plots
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, table, output_path):

        """
        This function ...
        :param table:
        :param output_path:
        :return:
        """

        # Call the setup of the base class
        super(Plotter, self).setup()

        # Inform the user
        log.info("Reading input data...")

        # Set the input table
        self.table = table

        # Set the path to the output directory
        self.config.output = output_path

    # -----------------------------------------------------------------

    @abstractmethod
    def prepare_data(self):

        """
        This function ...
        :return:
        """

        return

    # -----------------------------------------------------------------

    @abstractmethod
    def plot(self):

        """
        This function ...
        :return:
        """

        return

# -----------------------------------------------------------------
