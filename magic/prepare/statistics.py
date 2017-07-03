#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.prepare.statistics Contains the PreparationStatistics class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.composite import SimplePropertyComposite

# -----------------------------------------------------------------

class PreparationStatistics(SimplePropertyComposite):

    """
    This function ...
    """

    def __init__(self, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(PreparationStatistics, self).__init__()

        # Define properties
        self.add_property("convolution_filter", "filter", "convolution filter")
        self.add_property("rebinning_filter", "filter", "rebinning filter")
        self.add_property("not_rebinned", "string_list", "images that are not rebinned")
        self.add_property("not_convolved", "string_list", "images that are not convolved")

        # Set properties
        self.set_properties(kwargs)

# -----------------------------------------------------------------
