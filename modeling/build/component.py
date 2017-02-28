#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.build.component Contains the BuildComponent class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..component.galaxy import GalaxyModelingComponent
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

class BuildComponent(GalaxyModelingComponent):
    
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
        super(BuildComponent, self).__init__(config)

        # Paths
        self.model_stars_path = None
        self.model_dust_path = None

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(BuildComponent, self).setup()

        # Paths
        self.model_stars_path = fs.create_directory_in(self.model_path, "stars")
        self.model_dust_path = fs.create_directory_in(self.model_path, "dust")
 
# -----------------------------------------------------------------
