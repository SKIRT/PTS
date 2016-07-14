#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.decomposition.component Contains the DecompositionComponent class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..core.component import ModelingComponent
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

class DecompositionComponent(ModelingComponent):
    
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
        super(DecompositionComponent, self).__init__(config)

        # -- Attributes --

        # The path to the bulge, disk and model directories
        self.bulge_path = None
        self.disk_path = None
        self.model_path = None

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(DecompositionComponent, self).setup()

        # Determine the path to the bulge, disk and model directories
        self.bulge_path = fs.join(self.components_path, "bulge")
        self.disk_path = fs.join(self.components_path, "disk")
        self.model_path = fs.join(self.components_path, "model")

        # Create the bulge and disk directories
        fs.create_directories(self.bulge_path, self.disk_path, self.model_path)

# -----------------------------------------------------------------
