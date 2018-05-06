#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.energy.projected Contains the ProjectedEnergyAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ..component import AnalysisComponent
from ....core.tools import filesystem as fs
from ....core.basics.log import log

# -----------------------------------------------------------------

class ProjectedEnergyAnalyser(AnalysisComponent):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(ProjectedEnergyAnalyser, self).__init__(*args, **kwargs)

        # -- Attributes --

    # -----------------------------------------------------------------

# -----------------------------------------------------------------
