#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.prepare.preparer Contains the ImagePreparer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from ...core.basics.log import log

# -----------------------------------------------------------------

class ImagePreparer(Configurable):
    
    """
    This class ...
    """
    
    def __init__(self, *args, **kwargs):
        
        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ImagePreparer, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

   
# -----------------------------------------------------------------
