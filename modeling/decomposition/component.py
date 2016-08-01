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

        # The path to the components/parameters directory
        self.components_parameters_path = None

        # The path to the components/images directory
        self.components_images_path = None

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(DecompositionComponent, self).setup()

        # Set ...
        self.components_parameters_path = fs.create_directory_in(self.components_path, "parameters")

        # Set ...
        self.components_images_path = fs.create_directory_in(self.components_path, "images")

# -----------------------------------------------------------------
