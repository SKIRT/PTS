#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from pts.core.test.implementation import TestImplementation
from pts.core.basics.log import log
from pts.modeling.fitting.tables import GenerationsTable
from pts.modeling.fitting.explorer import GenerationInfo
from pts.core.basics.map import Map
from pts.core.tools import filesystem as fs
from pts.core.tools import formatting as fmt

# -----------------------------------------------------------------

description = "testing launching simulations on remote hosts"

# -----------------------------------------------------------------

class LaunchTest(TestImplementation):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(TablesTest, self).__init__(*args, **kwargs)

        # The created tables
        self.tables = dict()

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

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

    

# -----------------------------------------------------------------
