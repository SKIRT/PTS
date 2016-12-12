#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.python Contains the RemotePythonSession class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# -----------------------------------------------------------------

class SimulationQueue(object):

    """
    This class ...
    """

    def __init__(self, remote):

        """
        This function ...
        :param service:
        :return:
        """

        self.remote = remote
        
        self.pid = None

        self.pip_file_name = None
        
# -----------------------------------------------------------------
