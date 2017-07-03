#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.plot.fwhm Contains the FwhmPlotter class

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
import matplotlib.pyplot as plt

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable

# -----------------------------------------------------------------

class FwhmPlotter(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(FwhmPlotter, self).__init__(*args, **kwargs)

        # -- Attributes --

# -----------------------------------------------------------------
