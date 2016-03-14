#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.basics.trackrecord Contains the TrackRecord class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import defaultdict

# -----------------------------------------------------------------

class TrackRecord(object):

    """
    This class ...
    """

    def __init__(self, initial_stage):

        """
        The constructor ...
        :return:
        """

        self.stage = initial_stage
        self.snapshots = defaultdict(list)

    # -----------------------------------------------------------------

    def set_stage(self, stage):

        """
        This function ...
        """

        self.stage = stage

    # -----------------------------------------------------------------

    def append(self, source):

        """
        This function ...
        """

        self.snapshots[self.stage].append(source)

    # -----------------------------------------------------------------

    def plot(self):
        
        """
        This function ...
        """

        # Loop over all snapshots
        for title, stage in self.snapshots.items():

            # Plot ...
            for snapshot in stage: snapshot.plot(title=title)
    
# -----------------------------------------------------------------