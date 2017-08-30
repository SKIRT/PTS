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

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.tools import filesystem as fs
from ..analyse.database import load_database, get_runs
from ...core.tools.serialization import write_dict
from ...core.basics.log import log

# -----------------------------------------------------------------

class Plotter(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(Plotter, self).__init__(*args, **kwargs)

        # The database connection
        self.database = None

        # The plot
        #self.plot = None

        # The plot data
        self.data = None

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

    def set_title(self, title):

        """
        This function ...
        :param title: 
        :return: 
        """

        title = title.replace('_', '\_')
        plt.title(title)

    # -----------------------------------------------------------------

    def write_or_show(self, name=None):

        """
        This function ...
        :param name: 
        :return: 
        """

        if "output" in self.config:

            path = self.output_path_file(name + "." + self.config.format)
            plt.savefig(path)

        else: plt.show()

        plt.close()

    # -----------------------------------------------------------------

    def write_data(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Writing the data ...")

        # Determine the path
        path = self.output_path_file("data_" + self.class_name + ".dat")

        # Write the data
        write_dict(self.data, path)

    # -----------------------------------------------------------------

    def clear_data(self):

        """
        This function ...
        :return: 
        """

        self.data = None

    # -----------------------------------------------------------------

    @property
    def runs(self):

        """
        This function ...
        :return: 
        """

        if self.config.runs is None: return get_runs(self.database)
        else: return self.config.runs

# -----------------------------------------------------------------
