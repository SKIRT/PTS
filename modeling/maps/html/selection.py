#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.html.selection Contains the MapsSelectionPageGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import defaultdict

# Import the relevant PTS classes and modules
from ....core.basics.log import log
from ..component import MapMakingComponent

# -----------------------------------------------------------------

plots_name = "plots"
ncolumns = 2
background_color = "white"
key_color = "#4180d3"

# -----------------------------------------------------------------

page_width = 600

# -----------------------------------------------------------------

intervals_filename = "intervals.dat"

# -----------------------------------------------------------------

class MapsSelectionPageGenerator(MapMakingComponent):

    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(MapsSelectionPageGenerator, self).__init__(*args, **kwargs)

# -----------------------------------------------------------------
