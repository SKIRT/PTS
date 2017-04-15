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
from pts.core.basics.table import SmartTable
from pts.core.tools import tables

# -----------------------------------------------------------------

class ScoresTable(SmartTable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ScoresTable, self).__init__(*args, **kwargs)

        # Add column info
        self.column_info.append(("Individual name", str, None, "name of the individual"))
        self.column_info.append(("Score", float, None, "individual's score"))

    # -----------------------------------------------------------------

    @property
    def best_individual_name(self):

        """
        This function ...
        :return:
        """

        index = np.argmin(self["Score"])
        return self["Individual name"][index]

    # -----------------------------------------------------------------

    def score_for(self, individual_name):

        """
        This function ...
        :param individual_name: 
        :return: 
        """

        index = tables.find_index(self, individual_name, "Individual name")
        return self["Score"][index]

# -----------------------------------------------------------------