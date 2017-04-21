#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.test.implementation import TestImplementation
from pts.core.tools.logging import log
from pts.modeling.fitting.tables import GenerationsTable
from pts.modeling.fitting.explorer import GenerationInfo
from pts.core.basics.map import Map
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

description = "testing the tables"

# -----------------------------------------------------------------

class TablesTest(TestImplementation):

    """
    This class ...
    """

    def __init__(self, config=None, interactive=False):

        """
        This function ...
        :param config:
        :param interactive:
        """

        # Call the constructor of the base class
        super(TablesTest, self).__init__(config, interactive)

        # The created tables
        self.tables = dict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # Test
        self.create_tables()

        # Write
        self.write()

        # Load
        self.load_tables()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(TablesTest, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def create_tables(self):

        """
        This function ...
        :return: 
        """

        table = GenerationsTable(parameters=["a"], units={"a": "m"})

        info = GenerationInfo()

        table.add_entry(info, {"a": Map(min=None, max=None)})

        #print(tab)

        #print(type(table["Generation index"][0]))

        self.tables["generations"] = table

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Writing ...")

        # Write the tables
        self.write_tables()

    # -----------------------------------------------------------------

    def write_tables(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Writing the tables ...")

        # Loop over the tables
        for name in self.tables:

            # Determine the path
            path = fs.join(self.path, name + ".dat")

            # Write
            self.tables[name].saveto(path)

    # -----------------------------------------------------------------

    def load_tables(self):

        """
        This function ...
        :return: 
        """

        # Loop over the tables
        for name in self.tables:

            # Determine path
            path = fs.join(self.path, name + ".dat")

            # Load
            table = GenerationsTable.from_file(path)

            print(name)
            print(table)

# -----------------------------------------------------------------
