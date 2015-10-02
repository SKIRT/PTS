#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Import Python 3 functionality
from __future__ import (absolute_import, division, print_function)

# *****************************************************************

class TrackRecord(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        self.snapshots = []

    # *****************************************************************

    def plot(self):
        
        """
        This function ...
        """

        # Loop over all snapshots
        for snapshot in self.snapshots: snapshot.plot()
    
# *****************************************************************