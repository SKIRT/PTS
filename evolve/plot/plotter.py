#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.evolve.plot.plotter Contains the Plotter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import matplotlib.pyplot as plt
import sqlite3

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.basics.plot import Plot
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

class Plotter(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        """

        # Call the constructor of the base class
        super(Plotter, self).__init__(config)

        # The database connection
        self.database = None

        # The plot
        self.plot = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

# -----------------------------------------------------------------
