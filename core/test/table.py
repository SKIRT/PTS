#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.test.table Contains the PTSTestsTable class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import imp
import importlib
from collections import defaultdict, OrderedDict

# Import the relevant PTS classes and modules
from ..basics.log import log
from ..basics.configurable import Configurable
from ..tools import introspection
from ..tools import filesystem as fs
from ..basics.configuration import ConfigurationDefinition, InteractiveConfigurationSetter, DictConfigurationSetter, PassiveConfigurationSetter
from .imports import ImportsChecker
from .test import PTSTest
from ..tools import time
from ..tools import stringify
from ..remote.utils import DetachedCalculation
from ..basics.table import SmartTable

# -----------------------------------------------------------------

class PTSTestsTable(SmartTable):

    """
    This class ...
    """

    # Add column info
    _column_info = OrderedDict()
    _column_info["Name"] = (str, None, "name of the test")
    _column_info["Start time"] = (str, None, "timestamp for start of command")
    _column_info["End time"] = (str, None, "timestamp for end of command")

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(PTSTestsTable, self).__init__(*args, **kwargs)

        # Add column info
        self.add_all_column_info(self._column_info)

    # -----------------------------------------------------------------

    def add_test(self, test):

        """
        This function ...
        :param test:
        :return:
        """

        values = []
        self.add_row(values)

# -----------------------------------------------------------------

# Determine path, initialize if not present
tests_table_path = fs.join(introspection.pts_tests_dir, "tests.dat")
if not fs.is_file(tests_table_path):

    # Create the table and save it
    table = PTSTestsTable()
    table.saveto(tests_table_path)

# -----------------------------------------------------------------

def load_tests_table():

    """
    This function ...
    :return:
    """

    return PTSTestsTable.from_file(tests_table_path)

# -----------------------------------------------------------------
