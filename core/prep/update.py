#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.prep.update
#

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable

# -----------------------------------------------------------------

class SkirtUpdater(Configurable):
    
    """
    This class ...
    """
    
    def __init__(self):
        
        """
        The constructor ...
        """
        
        # Call the constructor of the base class
        super(SkirtUpdater, self).__init__()
        
    # -----------------------------------------------------------------
    
    def run(self):
        
        """
        This function ...
        """

        # 1. Call the setup function
        self.setup()

# -----------------------------------------------------------------
