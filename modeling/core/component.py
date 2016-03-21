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
from ...core.tools import inspection, filesystem

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

        # Attributes
        self.galaxy_name = None

        # Modeling directories
        self.data_path = None
        self.prep_path = None
        self.truncation_path = None
        self.phot_path = None
        self.maps_path = None
        self.components_path = None
        self.fit_path = None
        self.analysis_path = None

        # PTS directories
        self.kernels_path = None

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(ModelingComponent, self).setup()

        # -- Attributes --

        # Get the name of the galaxy (the name of the base directory)
        self.galaxy_name = os.path.basename(self.config.path)

        # Get the full paths to the necessary subdirectories
        self.data_path = filesystem.join(self.config.path, "data")
        self.prep_path = os.path.join(self.config.path, "prep")
        self.truncation_path = filesystem.join(self.config.path, "truncated")
        self.phot_path = filesystem.join(self.config.path, "phot")
        self.maps_path = os.path.join(self.config.path, "maps")
        self.components_path = os.path.join(self.config.path, "components")
        self.fit_path = os.path.join(self.config.path, "fit")
        self.analysis_path = os.path.join(self.config.path, "analysis")

        # Determine the path to the kernels user directory
        self.kernels_path = os.path.join(inspection.pts_user_dir, "kernels")

        # Create the prep path if it does not exist yet
        filesystem.create_directories([self.prep_path, self.truncation_path, self.maps_path, self.phot_path, self.maps_path, self.components_path, self.fit_path, self.analysis_path])

# -----------------------------------------------------------------
