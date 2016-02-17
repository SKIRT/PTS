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

# Import astronomical modules
from astropy.table import Table
from astropy.io import ascii

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable
from ..tools.logging import log

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
        super(Plotter, self).__init__(None, "core")

        # -- Attributes --

        # The table containing the input data
        self.table = None

        # The data structure to contain the input information in plottable format
        self.data = []

        # The path to the output directory
        self.output_path = None

    # -----------------------------------------------------------------

    @staticmethod
    @abstractproperty
    def fill_values():

        """
        This function ...
        :return:
        """

        return []

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

    def run(self, input, output_path):

        """
        This function should be called to invoke the plotting routine ...
        :param input:
        :param output_path:
        :return:
        """

        # 1. Call the setup function
        self.setup(input, output_path)

        # 2. Prepare the input data into plottable format
        self.prepare_data()

        # 3. Create the plots
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, input, output_path):

        """
        This function ...
        :return:
        """

        # Call the setup of the base class
        super(Plotter, self).setup()

        # Inform the user
        log.info("Reading input data...")

        # If the input is a Table object
        if isinstance(input, Table): self.table = input

        # If the input is a string
        elif isinstance(input, basestring): self.table = ascii.read(input, fill_values=self.fill_values)

        # Invalid input
        else: raise ValueError("Input must be either an Astropy Table object or a filename (e.g. " + self.default_input() + ")")

        # Set the path to the output directory
        self.output_path = output_path

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
        :param path:
        :return:
        """

        return

# -----------------------------------------------------------------
