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

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.basics.plot import Plot
from ...core.tools import filesystem as fs
from ..analyse.database import load_database

# -----------------------------------------------------------------

class Plotter(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        """

        # Call the constructor of the base class
        super(Plotter, self).__init__(config, interactive)

        # The database connection
        self.database = None

        # The plot
        self.plot = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs: 
        :return: 
        """

        # Call the setup function of the base class
        super(Plotter, self).setup(**kwargs)

        # Load the database
        if "database" in kwargs: self.database = kwargs.pop("database")
        elif self.config.database is not None: self.database = load_database(self.config.database)
        else: raise ValueError("Database not specified")

# -----------------------------------------------------------------
