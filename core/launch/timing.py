#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.launch.parallelization Contains the Parallelization class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..tools import tables
from ..tools import filesystem as fs

# -----------------------------------------------------------------

class TimingTable(object):

    """
    This class ...
    """

    def __init__(self, path):

        """
        This function ...
        """

        # Set the path of the timing table
        self.path = path

        # If the file does not exist yet, create it
        if not fs.is_file(self.path): self.initialize()

    # -----------------------------------------------------------------

    def initialize(self):

        """
        This function ...
        :return:
        """



    # -----------------------------------------------------------------

    def add_entry(self, name, timestamp, ):

# -----------------------------------------------------------------
