#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.build.builder Contains the ModelBuilder class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import BuildComponent

# -----------------------------------------------------------------

class ModelBuilder(BuildComponent):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(ModelBuilder, self).__init__(config)

        # The maps
        self.old_stars = None
        self.young_stars = None
        self.ionizing_stars = None
        self.dust = None

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(ModelBuilder, self).setup()

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """
    
        self.write_maps()

    # -----------------------------------------------------------------

    def write_maps(self):

        """
        This function ...
        :return:
        """

# -----------------------------------------------------------------
