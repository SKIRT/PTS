#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.component Contains the ModelingComponent class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.tools import inspection

# -----------------------------------------------------------------

class ModelingComponent(Configurable):
    
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
        super(ModelingComponent, self).__init__(config, "modeling")

    # -----------------------------------------------------------------

    def setup(self, path):

        """
        This function ...
        :return:
        """

        # -- Attributes --

        # Get the name of the galaxy (the name of the base directory)
        self.galaxy_name = os.path.basename(path)

        # Get the full path to the 'data', 'prep' and 'in' directories
        self.data_path = os.path.join(path, self.config.data_dir)
        self.prep_path = os.path.join(path, self.config.prep_dir)
        self.in_path = os.path.join(path, self.config.in_dir)

        # Determine the path to the kernels user directory
        self.kernels_path = os.path.join(inspection.pts_user_dir, "kernels")

# -----------------------------------------------------------------
